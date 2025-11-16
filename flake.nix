{
  description = "Rust dev flake";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    rust-overlay.url = "github:oxalica/rust-overlay";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = {
    self,
    nixpkgs,
    rust-overlay,
    flake-utils,
    ...
  }:
    flake-utils.lib.eachDefaultSystem (
      system: let
        overlays = [(import rust-overlay)];
        pkgs = import nixpkgs {
          inherit system overlays;
        };
      in {
        devShells.default = with pkgs;
          mkShell {
            buildInputs = [
              openssl
              pkg-config
              eza
              fd
              cargo-semver-checks
              rust-bin.stable.latest.default
            ];
          };
        packages.x86_64-linux.default = pkgs.rustPlatform.buildRustPackage {
          name = "";
          src = ./.;
          cargoLock.lockFile = ./Cargo.lock;
        };
      }
    );
}
