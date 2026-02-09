{
  description = "GAND Project - OCI image stable nix";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-23.11";
  };
  
  nixConfig = {
    extra-substituters = [ "https://cache.nixos.org" ];
    allow-import-from-derivation = true;
  };

  outputs = { self, nixpkgs }:
    let
      system = "x86_64-linux";
      pkgs = import nixpkgs { inherit system; };
      
      # Define the R environment with the required packages
      Renv = pkgs.rWrapper.override {
        packages = with pkgs.rPackages; [
          future Matrix Seurat dplyr tidyr 
          patchwork DESeq2 ggplot2 rmarkdown 
          hdf5r ggpubr lme4 emmeans 
          argparser jsonlite RColorBrewer ggtext 
          tinytex
        ];
      };
    in
    {
      packages.${system}.default = pkgs.dockerTools.streamLayeredImage {
        name = "gand_nix_image";
        tag = "v0.0.0";
        contents = [ 
          Renv 
          pkgs.bash 
          pkgs.coreutils 
          pkgs.pandoc
        ];
        config = {
          Cmd = [ "${Renv}/bin/R" ];
        };
      };
    };
}