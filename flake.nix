{
  description = "A very basic flake";

  outputs = { self, nixpkgs }: {

    env = nixpkgs.legacyPackages.x86_64-linux.bundlerEnv {
      name = "setser.github.io-bundler-env";
      ruby = nixpkgs.legacyPackages.x86_64-linux.ruby;
      gemfile  = ./Gemfile;
      lockfile = ./Gemfile.lock;
      gemset   = ./gemset.nix;
      groups   = [ "default" "jekyll_plugins" ];
    };

    defaultPackage.x86_64-linux = nixpkgs.legacyPackages.x86_64-linux.stdenv.mkDerivation {
      name = "setser.github.io";
      buildInputs = [ self.env ];
      src = ./.;
      installPhase = [ true ];
      shellHook = ''
        exec ${self.env}/bin/jekyll serve -P 8080 --watch
      '';
    };


    devShell.x86_64-linux = nixpkgs.legacyPackages.x86_64-linux.mkShell {
      name = "setser.github.io";
      buildInputs = [ self.env ];
    };
  };
}
