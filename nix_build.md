NOTE: This is experiemental and will be removed in the future unless 
a stable build process can be provided.

For highly reproducible environment builds, we also provide a `Nix` flake which will build a OCI/Docker image. Using Nix to build the docker image ensures that all package versions and underlying libraries will be as close to bit-for-bit reproducibl across machines and time. First, make sure that you have `Nix` installed. You can find it [here](https://nixos.org/download/). Note that `Nix` is not compatible with Windows but can work through `WSL`.

The resulting image with also use the `v0.0.1` nomenclature but the image will be called `gand_nix_image.tar` to clarify that this is strictly built image less prone to config drift. 

# Nix on Linux
Once `Nix` is installed, run the following:

```
cd containers
nix build .#default
./result > GAND_nix_image.tar
```


# Nix OCI from mac

==NOTE: this section is not easy to make work and at the moments fails more often than not==

To build the image with `Nix` on Mac using Linux as target:

Run this once on the Mac:

1. `sudo nano /etc/nix/nix.conf`

2. Add this line: `trusted-users = root @admin` (This trusts all Mac admins).

3. Restart Nix: `sudo launchctl kickstart -k system/org.nixos.nix-daemon`

Open terminal to run the nix Deamon:

```
cd containers
nix run nixpkgs#darwin.linux-builder-x86_64
```
Then in your another terminal run the build process. You will be prompted to allow config settings. Type `y` to accept the settings.

```
nix build .#packages.x86_64-linux.default \
  --option builders "ssh://builder x86_64-linux - 4 1 kvm,nixos-test,benchmark,big-parallel" \
  --max-jobs 0

./result > GAND_image.tar
```

To exit the `Nix` builder, run in that terminal:

```
shutdown now
```

