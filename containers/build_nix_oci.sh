#!/bin/bash
set -e

# 1. Check if the Linux builder is reachable
if ! ssh -q builder exit; then
  echo "❌ Error: Linux builder is not running."
  echo "Please run: nix run nixpkgs#darwin.linux-builder-x86_64 in another tab."
  exit 1
fi

echo "🚀 Building Linux OCI image..."

# 2. Run the build with the correct flags
# We use --max-jobs 0 to force the Mac to send the work to the Linux builder
nix build .#packages.x86_64-linux.default \
  --option builders "ssh://builder x86_64-linux" \
  --max-jobs 0

# 3. Export the result to a tarball
echo "📦 Exporting to GAND_image.tar..."
./result > GAND_image.tar

echo "✅ Done! You can now transfer GAND_image.tar to the HPC."