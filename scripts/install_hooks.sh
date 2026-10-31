#!/bin/sh
# Install a git pre-commit hook that runs scripts/leak_check.py on the staged files and blocks the commit on hits.
#   scripts/install_hooks.sh            install (an existing foreign hook is kept as pre-commit.local and still run)
#   scripts/install_hooks.sh --uninstall
# The hook reads the private data folder from $DNA_DATA (see leak_check.py). Bypass once with: git commit --no-verify
set -eu

ROOT=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
HOOKS=$(git -C "$ROOT" rev-parse --git-path hooks)
case "$HOOKS" in /*) ;; *) HOOKS="$ROOT/$HOOKS" ;; esac
HOOK="$HOOKS/pre-commit"
TAG="# installed by scripts/install_hooks.sh (leak_check)"

if [ "${1:-}" = "--uninstall" ]; then
  if [ -f "$HOOK" ] && grep -qF "$TAG" "$HOOK"; then
    rm "$HOOK"
    [ -f "$HOOKS/pre-commit.local" ] && mv "$HOOKS/pre-commit.local" "$HOOK"
    echo "Removed the leak_check pre-commit hook."
  else
    echo "No leak_check pre-commit hook installed."
  fi
  exit 0
fi

mkdir -p "$HOOKS"
if [ -f "$HOOK" ] && ! grep -qF "$TAG" "$HOOK"; then
  mv "$HOOK" "$HOOKS/pre-commit.local"
  echo "Kept the existing pre-commit hook as pre-commit.local (it still runs after the leak check)."
fi

cat > "$HOOK" <<EOF
#!/bin/sh
$TAG
root=\$(git rev-parse --show-toplevel)
py=\${PYTHON:-python3}
if ! "\$py" "\$root/scripts/leak_check.py" --staged --repo "\$root"; then
  echo "pre-commit: commit blocked by leak_check (see above). Bypass only if sure: git commit --no-verify" >&2
  exit 1
fi
if [ -x "\$(dirname "\$0")/pre-commit.local" ]; then
  exec "\$(dirname "\$0")/pre-commit.local" "\$@"
fi
EOF
chmod +x "$HOOK"
echo "Installed pre-commit hook: $HOOK"
