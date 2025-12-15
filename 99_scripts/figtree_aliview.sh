#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# Config (override via env vars)
# -----------------------------
FIGTREE_VERSION="${FIGTREE_VERSION:-1.4.4}"
ALIVIEW_VERSION="${ALIVIEW_VERSION:-1.30}"

BASE_DIR="${BASE_DIR:-$HOME/.local/opt}"
BIN_DIR="${BIN_DIR:-$HOME/.local/bin}"

# -----------------------------
# Helpers
# -----------------------------
have_cmd() { command -v "$1" >/dev/null 2>&1; }

download() {
  local url="$1"
  local out="$2"
  if have_cmd curl; then
    curl -L --fail -o "$out" "$url"
  elif have_cmd wget; then
    wget -O "$out" "$url"
  else
    echo "Error: need curl or wget." >&2
    exit 1
  fi
}

msg() { echo "[INFO] $*"; }

TMP_DIR="$(mktemp -d)"
cleanup() { rm -rf "$TMP_DIR"; }
trap cleanup EXIT

mkdir -p "$BASE_DIR" "$BIN_DIR"

# -----------------------------
# Install FigTree
# -----------------------------
FIGTREE_URL="https://github.com/rambaut/figtree/releases/download/v${FIGTREE_VERSION}/FigTree_v${FIGTREE_VERSION}.tgz"
FIGTREE_ARCHIVE="$TMP_DIR/figtree.tgz"
FIGTREE_DIR="$BASE_DIR/figtree"

msg "Downloading FigTree ${FIGTREE_VERSION}..."
download "$FIGTREE_URL" "$FIGTREE_ARCHIVE"

msg "Installing FigTree into ${FIGTREE_DIR}..."
rm -rf "$FIGTREE_DIR"
mkdir -p "$FIGTREE_DIR"
tar -xzf "$FIGTREE_ARCHIVE" -C "$FIGTREE_DIR"

# Create FigTree launcher
cat > "$BIN_DIR/figtree" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

FIGTREE_DIR="${FIGTREE_DIR:-$HOME/.local/opt/figtree}"

# Prefer a packaged launcher if present
if [ -x "$FIGTREE_DIR/FigTree" ]; then
  exec "$FIGTREE_DIR/FigTree" "$@"
elif [ -x "$FIGTREE_DIR/figtree" ]; then
  exec "$FIGTREE_DIR/figtree" "$@"
fi

# Otherwise fall back to running the jar
JAR=""
if [ -f "$FIGTREE_DIR/lib/figtree.jar" ]; then
  JAR="$FIGTREE_DIR/lib/figtree.jar"
else
  JAR="$(find "$FIGTREE_DIR" -maxdepth 4 -type f -iname '*figtree*.jar' 2>/dev/null | head -n1 || true)"
fi

if [ -z "$JAR" ]; then
  echo "Could not find FigTree launcher or jar in $FIGTREE_DIR" >&2
  exit 1
fi

if ! command -v java >/dev/null 2>&1; then
  echo "Java is required to run FigTree but was not found in PATH." >&2
  echo "Ask your administrator to provide a Java runtime or load a module if available." >&2
  exit 1
fi

exec java -jar "$JAR" "$@"
EOF
chmod +x "$BIN_DIR/figtree"

# -----------------------------
# Install AliView
# -----------------------------
ALIVIEW_URL="https://ormbunkar.se/aliview/downloads/linux/linux-version-${ALIVIEW_VERSION}/aliview.tgz"
ALIVIEW_ARCHIVE="$TMP_DIR/aliview.tgz"
ALIVIEW_DIR="$BASE_DIR/aliview"

msg "Downloading AliView ${ALIVIEW_VERSION}..."
download "$ALIVIEW_URL" "$ALIVIEW_ARCHIVE"

msg "Installing AliView into ${ALIVIEW_DIR}..."
rm -rf "$ALIVIEW_DIR"
mkdir -p "$ALIVIEW_DIR"
tar -xzf "$ALIVIEW_ARCHIVE" -C "$ALIVIEW_DIR"

# Ensure bundled script is executable if present
if [ -f "$ALIVIEW_DIR/aliview" ]; then
  chmod +x "$ALIVIEW_DIR/aliview"
fi

# Create AliView launcher
cat > "$BIN_DIR/aliview" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

ALIVIEW_DIR="${ALIVIEW_DIR:-$HOME/.local/opt/aliview}"

# Prefer the provided start script
if [ -x "$ALIVIEW_DIR/aliview" ]; then
  exec "$ALIVIEW_DIR/aliview" "$@"
fi

# Fallback to jar
JAR=""
if [ -f "$ALIVIEW_DIR/aliview.jar" ]; then
  JAR="$ALIVIEW_DIR/aliview.jar"
else
  JAR="$(find "$ALIVIEW_DIR" -maxdepth 3 -type f -iname 'aliview*.jar' 2>/dev/null | head -n1 || true)"
fi

if [ -z "$JAR" ]; then
  echo "Could not find AliView script or jar in $ALIVIEW_DIR" >&2
  exit 1
fi

if ! command -v java >/dev/null 2>&1; then
  echo "Java is required to run AliView but was not found in PATH." >&2
  exit 1
fi

exec java -jar "$JAR" "$@"
EOF
chmod +x "$BIN_DIR/aliview"

# -----------------------------
# PATH hint
# -----------------------------
msg "Done."
msg "Launchers created:"
msg "  $BIN_DIR/figtree"
msg "  $BIN_DIR/aliview"
msg ""
msg "If needed, enable them in this shell:"
msg "  export PATH=\"$HOME/.local/bin:\$PATH\""
msg ""
msg "To make it permanent:"
msg "  echo 'export PATH=\"$HOME/.local/bin:\$PATH\"' >> ~/.bashrc"
 
