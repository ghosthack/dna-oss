"""Tests for the steps the offline demo run skips: polygenic scores, 1000 Genomes ancestry, Neanderthal and
fine-scale (AADR) ancestry. They use small hand-built inputs, so no network, bcftools or real data is needed.

Run from the repo root:  python3 -m unittest discover -s tests
ancestry_ho needs numpy (requirements.txt); its tests are skipped without it.
"""
import gzip
import math
import os
import sys
import tempfile
import unittest

os.environ.setdefault("DNA_DATA", tempfile.mkdtemp(prefix="dna-oss-test-"))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))

import ancestry  # noqa: E402
import neanderthal  # noqa: E402
import pgs  # noqa: E402

try:
    import ancestry_ho  # noqa: E402
    import ancestry_ho_fetch  # noqa: E402
except ImportError:  # numpy missing
    ancestry_ho = None


class ScoreAllTest(unittest.TestCase):
    """pgs.score_all on a hand-built 1000 Genomes cache: 4 European and 2 American reference samples."""

    SAMPLES = ["E1", "E2", "E3", "E4", "A1", "A2"]

    def setUp(self):
        self.out = tempfile.mkdtemp()
        ref = os.path.join(self.out, "ref")
        os.makedirs(ref)
        with open(os.path.join(ref, "samples.txt"), "w") as fh:
            fh.write("\n".join(self.SAMPLES) + "\n")
        with open(os.path.join(ref, "ped.txt"), "w") as fh:
            fh.write("FamilyID SampleID FatherID MotherID Sex Population Superpopulation\n")
            for s in self.SAMPLES:
                pop, sup = ("IBS", "EUR") if s.startswith("E") else ("MXL", "AMR")
                fh.write(f"{s} {s} 0 0 1 {pop} {sup}\n")
        # chrom pos ref alt, then one GT per sample (same order as samples.txt)
        rows = [
            ("chr1", 100, "A", "G", ["0|0", "0|1", "1|1", "0|1", "0|0", "1|1"]),    # plain SNP
            ("chr1", 200, "A", "T", ["0|0", "0|0", "0|0", "0|1", "0|0", "0|0"]),    # A/T, T rare in EUR (12.5%)
            ("chr1", 300, "C", "G", ["1|1", "1|1", "1|1", "0|1", "1|1", "1|1"]),    # C/G, G common in EUR (87.5%)
            ("chr1", 400, "A", "T", ["0|1"] * 6),                                   # A/T, no frequency in score
            ("chr1", 500, "C", "T", ["0|1", ".", "0|0", "0|0", "0|0", "0|0"]),      # missing reference genotype
        ]
        with gzip.open(os.path.join(ref, "chr1.tsv.gz"), "wt") as fh:
            for c, p, r, a, gts in rows:
                fh.write("\t".join([c, str(p), r, a] + gts) + "\n")
        self._saved = (pgs.OUT, list(pgs.GROUPS), pgs.LOCAL)
        pgs.OUT = self.out

    def tearDown(self):
        pgs.OUT, groups, pgs.LOCAL = self._saved
        pgs.GROUPS[:] = groups

    def score(self, variants):
        scores = {"PGS1": {"trait": "Test trait", "variants": variants, "published": len(variants)}}
        # The person: alt-allele dosages from imputation, {(chrom, pos): [(ref, alt, alt_dosage)]}
        person = {("chr1", 100): [("A", "G", 1.0)], ("chr1", 200): [("A", "T", 2.0)],
                  ("chr1", 300): [("C", "G", 2.0)], ("chr1", 400): [("A", "T", 1.0)],
                  ("chr1", 500): [("C", "T", 1.0)]}
        return pgs.score_all(scores, {"me": person})["me"]["PGS1"]

    def test_strand_handling_and_missing_reference(self):
        r = self.score([
            ("chr1", 100, "G", "A", 1.0, None),      # effect = ALT, dosage 1        -> 1
            ("chr1", 200, "T", "A", 10.0, 0.2),      # T reported 20%, EUR 12.5%: same strand -> 2 copies -> 20
            ("chr1", 300, "C", "G", 100.0, 0.8),     # C reported 80%, but EUR C is 12.5%: other strand, so the
                                                     # effect allele is really G (ALT), dosage 2 -> 200
            ("chr1", 400, "T", "A", 1000.0, None),   # palindromic, no frequency: can't tell the strand -> dropped
            ("chr1", 500, "T", "C", 10000.0, None),  # a reference genotype is missing -> skipped for everyone
        ])
        self.assertEqual(r["score"], 221.0)
        self.assertEqual(r["variants_used"], 3)
        self.assertEqual(r["variants_strand_ambiguous"], 1)

    def test_reference_distribution_uses_the_same_alleles(self):
        r = self.score([("chr1", 100, "G", "A", 1.0, None)])
        # EUR G counts are 0, 1, 2, 1 (mean 1); the person has 1 copy: one sample strictly below, z = 0
        self.assertEqual(r["groups"]["EUR"]["n"], 4)
        self.assertEqual(r["groups"]["ALL"]["n"], 6)
        self.assertEqual(r["groups"]["EUR"]["percentile"], 25.0)
        self.assertEqual(r["groups"]["EUR"]["z"], 0.0)

    def test_score_with_no_reference_spread_is_left_out(self):
        scores = {"PGS1": {"trait": "T", "variants": [("chr1", 400, "A", "T", 1.0, 0.1)], "published": 1}}
        # the only variant is palindromic with an unusable frequency, so nothing is scored
        self.assertNotIn("PGS1", pgs.score_all(scores, {"me": {("chr1", 400): [("A", "T", 1.0)]}})["me"])

    def test_unknown_local_population_is_dropped(self):
        pgs.LOCAL = "XXX"
        pgs.GROUPS[:] = [("XXX", "Typo")] + [g for g in pgs.GROUPS if g[0] != "XXX"]
        r = self.score([("chr1", 100, "G", "A", 1.0, None)])
        self.assertNotIn("XXX", r["groups"])
        self.assertNotIn("XXX", [c for c, _ in pgs.GROUPS])


class LoadScoresTest(unittest.TestCase):
    def test_reads_effect_allele_frequency(self):
        d = tempfile.mkdtemp()
        with gzip.open(os.path.join(d, "PGS9_hmPOS_GRCh38.txt.gz"), "wt") as fh:
            fh.write("#pgs_id=PGS9\n#trait_reported=Test\n#variants_number=3\n")
            fh.write("rsID\teffect_allele\tother_allele\teffect_weight\tallelefrequency_effect\thm_chr\thm_pos\n")
            fh.write("rs1\ta\tg\t0.5\t0.3\t1\t100\n")
            fh.write("rs2\tT\tA\t0.1\tNA\t2\t200\n")
            fh.write("rs3\tC\tT\t0.1\t0.1\tX\t300\n")  # chrX is not scored
        saved = pgs.PGS_DIR
        pgs.PGS_DIR = d
        try:
            s = pgs.load_scores()["PGS9"]
        finally:
            pgs.PGS_DIR = saved
        self.assertEqual(s["variants"], [("chr1", 100, "A", "G", 0.5, 0.3), ("chr2", 200, "T", "A", 0.1, None)])
        self.assertEqual(s["published"], 3)


class JackknifeTest(unittest.TestCase):
    def test_equal_groups_match_the_classic_jackknife(self):
        reps = [0.30, 0.34, 0.28, 0.33]
        m = sum(reps) / 4
        classic = math.sqrt(3 / 4 * sum((r - m) ** 2 for r in reps))
        self.assertAlmostEqual(ancestry.weighted_jackknife_se(reps, [10, 10, 10, 10], 0.31), classic)

    def test_unequal_groups_by_hand(self):
        # sizes 1 and 3 (n = 4): h = 4 and 4/3; est = 2*full - (0.75*r1 + 0.25*r2)
        full, r1, r2 = 0.5, 0.4, 0.6
        est = 2 * full - (0.75 * r1 + 0.25 * r2)
        var = ((4 * full - 3 * r1 - est) ** 2 / 3 + ((4 / 3) * full - (1 / 3) * r2 - est) ** 2 / (1 / 3)) / 2
        self.assertAlmostEqual(ancestry.weighted_jackknife_se([r1, r2], [1, 3], full), math.sqrt(var))

    def test_needs_two_groups(self):
        self.assertIsNone(ancestry.weighted_jackknife_se([0.3], [10], 0.3))

    def test_missing_genotype_keeps_samples_aligned(self):
        fix = lambda s: ancestry.GT_TOKEN.sub(lambda m: m.group(1) + m.group(2) if m.group(1) else "..", s)
        self.assertEqual(fix("0|1.1/1"), "01..11")


class NeanderthalTest(unittest.TestCase):
    def test_reference_dosages(self):
        self.assertEqual(neanderthal.ref_dosages(["0|1", "1|1", "0|0"], [0, 1, 2], True), [1, 2, 0])
        self.assertEqual(neanderthal.ref_dosages(["0|1", "1|1", "0|0"], [0, 1, 2], False), [1, 0, 2])

    def test_missing_reference_genotype_skips_the_tag(self):
        self.assertIsNone(neanderthal.ref_dosages(["0|1", "."], [0, 1], True))
        self.assertEqual(neanderthal.ref_dosages(["0|1", "."], [0], True), [1])  # only kept samples matter


@unittest.skipIf(ancestry_ho is None, "numpy not installed")
class FineScaleStrandTest(unittest.TestCase):
    """AADR alleles A/C on the opposite strand to the chips: ancestry_ho_fetch writes snps.tsv as T/G."""

    A1, A2 = "A".translate(ancestry_ho_fetch.COMP), "C".translate(ancestry_ho_fetch.COMP)  # "T", "G"

    def test_late_added_person_matches_the_cached_family(self):
        for chip in ("TT", "GT", "GG"):
            cached = ancestry_ho_fetch.dosage(chip, self.A1, False)  # how the fetch step stores family dosages
            self.assertEqual(ancestry_ho.chip_dosage(chip, self.A1, self.A2), cached, chip)

    def test_flipped_snp_is_not_dropped(self):
        self.assertEqual(ancestry_ho.chip_dosage("GT", self.A1, self.A2), 1)
        self.assertEqual(ancestry_ho.chip_dosage("AC", self.A1, self.A2), 3)  # AADR-strand letters don't match

    def test_grandpa_allele(self):
        self.assertEqual(ancestry_ho.haplotype_code("T", self.A1, self.A2), 1)
        self.assertEqual(ancestry_ho.haplotype_code("G", self.A1, self.A2), 0)
        self.assertEqual(ancestry_ho.haplotype_code("A", self.A1, self.A2), 3)


if __name__ == "__main__":
    unittest.main()
