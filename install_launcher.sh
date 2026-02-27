#!/usr/bin/env bash
# install_launcher.sh — One-time setup for the scRNA-seq desktop shortcut.
#
# Run this script ONCE from a terminal:
#
#   bash /home/yuta/Git/10Xapp/install_launcher.sh
#
# After this, you can double-click "scRNA-seq Analysis" on your Desktop
# to launch the app — no terminal or RStudio needed.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DESKTOP_DIR="$HOME/Desktop"
DESKTOP_FILE="$SCRIPT_DIR/10Xapp.desktop"
SHORTCUT="$DESKTOP_DIR/10Xapp.desktop"

echo ""
echo "=== scRNA-seq Analysis App — Desktop Setup ==="
echo ""

# ---- Validate files exist --------------------------------------------------

if [[ ! -f "$DESKTOP_FILE" ]]; then
  echo "ERROR: 10Xapp.desktop not found in $SCRIPT_DIR"
  echo "  Make sure you are running this script from inside the app folder."
  exit 1
fi

if [[ ! -f "$SCRIPT_DIR/launch.sh" ]]; then
  echo "ERROR: launch.sh not found in $SCRIPT_DIR"
  exit 1
fi

# ---- Make scripts executable -----------------------------------------------

chmod +x "$SCRIPT_DIR/launch.sh"
echo "[1/4] Made launch.sh executable."

# ---- Create Desktop directory if needed ------------------------------------

mkdir -p "$DESKTOP_DIR"

# ---- Generate .desktop file with the correct absolute Exec path -------------
# The repo copy of 10Xapp.desktop has a placeholder Exec path.
# sed rewrites it to the actual install location so the shortcut works
# regardless of where the project folder lives on this machine.

sed "s|^Exec=.*|Exec=$SCRIPT_DIR/launch.sh|" "$DESKTOP_FILE" > "$SHORTCUT"
chmod +x "$SHORTCUT"
echo "[2/4] Copied shortcut to $SHORTCUT (Exec set to $SCRIPT_DIR/launch.sh)"

# ---- Mark the shortcut as trusted (required on Ubuntu/GNOME) ---------------

if command -v gio &>/dev/null; then
  gio set "$SHORTCUT" metadata::trusted true
  echo "[3/4] Marked shortcut as trusted (gio)."
elif command -v dbus-launch &>/dev/null; then
  # Fallback for some desktop environments
  dbus-launch gio set "$SHORTCUT" metadata::trusted true 2>/dev/null \
    && echo "[3/4] Marked shortcut as trusted (dbus-launch gio)." \
    || echo "[3/4] Could not mark as trusted automatically — right-click the icon and choose 'Allow Launching'."
else
  echo "[3/4] Could not mark as trusted automatically."
  echo "      Right-click the icon on your Desktop and choose 'Allow Launching'."
fi

# ---- Done ------------------------------------------------------------------

echo "[4/4] Setup complete!"
echo ""
echo "You can now double-click 'scRNA-seq Analysis' on your Desktop to start the app."
echo "The first launch will install required R packages (10–30 min, one time only)."
echo ""
echo "Tip: You can also search for 'scRNA-seq' in your application launcher."
echo ""
