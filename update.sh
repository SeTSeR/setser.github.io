#!/usr/bin/env sh

nix shell nixpkgs#bundler nixpkgs#bundix nixpkgs#gcc nixpkgs#pkgconfig nixpkgs#gnumake nixpkgs#binutils nixpkgs#file -c bundler update
nix shell nixpkgs#bundler nixpkgs#bundix nixpkgs#gcc nixpkgs#pkgconfig nixpkgs#gnumake nixpkgs#binutils nixpkgs#file -c bundler lock
nix shell nixpkgs#bundler nixpkgs#bundix nixpkgs#gcc nixpkgs#pkgconfig nixpkgs#gnumake nixpkgs#binutils nixpkgs#file -c bundler package --no-install --path vendor
nix shell nixpkgs#bundler nixpkgs#bundix nixpkgs#gcc nixpkgs#pkgconfig nixpkgs#gnumake nixpkgs#binutils nixpkgs#file -c bundix
rm -fr vendor
