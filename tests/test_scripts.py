"""Unit tests for the analysis scripts (standard library only, no data or network needed).

Run from the repo root:  python3 -m unittest discover -s tests
"""
import gzip
import os
import sys
import tempfile
import unittest

# Point the scripts at an empty data folder before importing them (they read DNA_DATA at import time).
os.environ["DNA_DATA"] = tempfile.mkdtemp(prefix="dna-oss-test-")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))

import ancestry  # noqa: E402
import dnalib  # noqa: E402
import doctor  # noqa: E402
import hereditary_markers  # noqa: E402
import pgs  # noqa: E402
import pgx  # noqa: E402
import report_archive  # noqa: E402
import sibling_diffs  # noqa: E402
import wellness  # noqa: E402


def chip(**genotypes):
    """{rsid: (chrom, pos, genotype)} as dnalib.load_person returns it."""
    return {rs: ("1", 1, gt) for rs, gt in genotypes.items()}


class PgxTest(unittest.TestCase):
    def test_untested_single_snp_gene_is_not_normal(self):
        for gene in ("SLCO1B1", "VKORC1", "ABCG2", "CYP2B6", "CYP4F2", "HLA-B*57:01"):
            dip, pheno, _ = pgx.call(gene, {}, male=False)
            self.assertEqual(dip, "not covered", gene)
            self.assertTrue(pheno.startswith("Not determinable"), gene)

    def test_partial_panel_names_untested_alleles(self):
        _, pheno, _ = pgx.call("CYP2C19", chip(rs4244285="GG"), male=False)
        self.assertTrue(pheno.startswith("Normal metabolizer"))
        self.assertIn("*3, *17 not tested", pheno)

    def test_full_panel_has_no_caveat(self):
        _, pheno, _ = pgx.call("SLCO1B1", chip(rs4149056="CC"), male=False)
        self.assertEqual(pheno, "Poor function")

    def test_g6pd_a_minus_without_376_probe(self):
        _, pheno, _ = pgx.call("G6PD", {"rs1050828": ("X", 1, "T")}, male=True)
        self.assertTrue(pheno.startswith("Deficient"))

    def test_g6pd_376a_background_is_not_a_minus(self):
        _, pheno, _ = pgx.call("G6PD", {"rs1050828": ("X", 1, "T"), "rs1050829": ("X", 1, "T")}, male=True)
        self.assertTrue(pheno.startswith("Normal"))

    def test_g6pd_compound_heterozygous_female(self):
        data = chip(rs1050828="CT", rs1050829="CT", rs5030868="AG")
        _, pheno, _ = pgx.call("G6PD", data, male=False)
        self.assertTrue(pheno.startswith("Deficient"))

    def test_cyp2d6_star4_absorbs_star10(self):
        dip, _, _ = pgx.call("CYP2D6", chip(rs3892097="CT", rs1065852="AG", rs5030655="II", rs28371725="CC"),
                             male=False)
        self.assertEqual(dip, "*1/*4")


class WellnessTest(unittest.TestCase):
    def test_factor_v_leiden_reported_without_prothrombin_probe(self):
        g = {"rs6025": "CT", "rs1799963": None}
        self.assertTrue(wellness.report_ok("Health", "Inherited clotting", g))
        res = wellness.clotting(g, "F")
        self.assertIn("Factor V Leiden ×1", res)
        self.assertIn("prothrombin G20210A not tested", res)

    def test_negative_clotting_with_missing_probe_says_so(self):
        res = wellness.clotting({"rs6025": "CC", "rs1799963": None})
        self.assertTrue(res.startswith("No Factor V Leiden"))
        self.assertIn("not tested", res)

    def test_hla_b27_without_cross_check(self):
        self.assertTrue(wellness.report_ok("Health", "HLA-B27 (spine arthritis)",
                                           {"rs4349859": "AG", "rs116488202": None}))


class DoctorTest(unittest.TestCase):
    def test_apoe_not_determined(self):
        line = doctor.apoe_line({"apoe": "", "source": "not determinable"})
        self.assertIn("not determined", line)
        self.assertNotIn("average risk", line)

    def test_apoe_e3e3(self):
        self.assertIn("average risk", doctor.apoe_line({"apoe": "ε3/ε3", "source": "measured"}))


class EscapingTest(unittest.TestCase):
    def test_sibling_table_headers_are_escaped(self):
        out = sibling_diffs.table(["<img src=x onerror=alert(1)>"], [])
        self.assertNotIn("<img", out)
        self.assertIn("&lt;img", out)


class PgsTest(unittest.TestCase):
    def test_palindromic(self):
        self.assertTrue(pgs.palindromic(("A", "T", 0), "A", "T"))
        self.assertTrue(pgs.palindromic(("C", "G", 0), "G", ""))
        self.assertFalse(pgs.palindromic(("A", "G", 0), "A", "G"))

    def test_strand_check(self):
        self.assertEqual(pgs.strand_check(0.2, 0.25), "keep")
        self.assertEqual(pgs.strand_check(0.2, 0.75), "flip")
        self.assertIsNone(pgs.strand_check(0.45, 0.2))  # too close to 50% to tell
        self.assertIsNone(pgs.strand_check(None, 0.2))

    def test_alt_count_missing(self):
        self.assertEqual(pgs.alt_count("0|1"), 1)
        self.assertEqual(pgs.alt_count("1/1"), 2)
        self.assertIsNone(pgs.alt_count("."))
        self.assertIsNone(pgs.alt_count(".|."))


class CacheTest(unittest.TestCase):
    def test_append_gzip_keeps_members(self):
        path = os.path.join(tempfile.mkdtemp(), "c.tsv.gz")
        dnalib.append_gzip(path, "a\n")
        dnalib.append_gzip(path, "b\n")
        with gzip.open(path, "rt") as fh:
            self.assertEqual(fh.read(), "a\nb\n")
        self.assertFalse(os.path.exists(path + ".tmp"))


class ParsingTest(unittest.TestCase):
    def test_reference_genotype_tokens_keep_alignment(self):
        sub = ancestry.GT_TOKEN.sub
        fix = lambda s: sub(lambda m: m.group(1) + m.group(2) if m.group(1) else "..", s)
        self.assertEqual(fix("0|1.1|1./."), "01..11..")

    def test_report_file_without_date(self):
        self.assertEqual(report_archive.capture_date("23andMe_reports_summary_alex.pdf"), "")
        self.assertEqual(report_archive.capture_date("23andMe_reports_summary_2024-05-01.pdf"), "2024-05-01")


class ChipMergeTest(unittest.TestCase):
    def test_discordant_snps_dropped_and_strand_flips_counted(self):
        a = {f"rs{i}": ("1", i, "AG") for i in range(100)}
        b = dict(a)
        b["rs1"] = ("1", 1, "CT")  # the same call on the opposite strand
        merged, stats = dnalib.merge(a, b)
        self.assertNotIn("rs1", merged)
        self.assertEqual(stats["discordant"], 1)
        self.assertEqual(stats["strand_complement"], 1)

    def test_hemizygous_x_coding_is_not_discordant(self):
        _, stats = dnalib.merge({"rs1": ("X", 1, "A")}, {"rs1": ("X", 1, "AA")})
        self.assertEqual(stats["discordant"], 0)


class HereditaryMarkerTest(unittest.TestCase):
    def test_indel_probe_counts(self):
        self.assertEqual(hereditary_markers.copies("DI", "D", "i4000377", "p"), "1")
        self.assertEqual(hereditary_markers.copies("II", "D", "i4000377", "p"), "0")

    def test_unexpected_encoding_is_untested_not_zero(self):
        self.assertEqual(hereditary_markers.copies("AA", "D", "i4000377", "p"), "")

    def test_no_call(self):
        self.assertEqual(hereditary_markers.copies(None, "C", "rs34612342", "p"), "")


class DisclaimerTest(unittest.TestCase):
    def test_every_page_generator_carries_the_notice(self):
        scripts = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts")
        for name in sorted(os.listdir(scripts)):
            if not name.endswith(".py"):
                continue
            src = open(os.path.join(scripts, name), encoding="utf-8").read()
            if "<!doctype html>" in src.lower() or "<html" in src.lower():
                self.assertTrue("DISCLAIMER_HTML" in src or "DISCLAIMER_LEAD" in src, name)


if __name__ == "__main__":
    unittest.main()
