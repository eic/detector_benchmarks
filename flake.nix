{
  inputs.epic-nix.url = "github:veprbl/epic-nix";

  inputs.nixpkgs.follows = "epic-nix/nixpkgs";

  outputs = { self, epic-nix, nixpkgs }:
    let

      inherit (nixpkgs) lib;
      supportedSystems = [ "x86_64-linux" "x86_64-darwin" ];

    in
    {

      devShells = lib.genAttrs supportedSystems
        (system:
          with import nixpkgs {
            inherit system;
            overlays = [ epic-nix.overlays.default ];
          };
          {
            default =
              mkShell {
                buildInputs =
                  (builtins.attrValues epic-nix.packages.${system})
                  ++ (with epic-nix.packages.${system}; [
                    geant4.data.G4EMLOW
                    geant4.data.G4ENSDFSTATE
                    geant4.data.G4ENSDFSTATE
                    geant4.data.G4PARTICLEXS
                    geant4.data.G4PhotonEvaporation
                  ])
                  ++ [
                    nlohmann_json
                    snakemake
                    python3
                    python3Packages.awkward
                    python3Packages.boto3
                    python3Packages.bokeh
                    python3Packages.dask
                    python3Packages.dask-awkward
                    python3Packages.dask-histogram
                    python3Packages.distributed
                    python3Packages.matplotlib
                    python3Packages.notebook
                    python3Packages.numpy
                    python3Packages.pyhepmc
                    python3Packages.scipy
                    python3Packages.setuptools
                    python3Packages.uproot
                  ];
                shellHook = ''
                  unset JUPYTER_CONFIG_DIR JUPYTER_PATH
                  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${lib.makeLibraryPath [ fmt ]}
                  export S3_ACCESS_KEY="eicS3read"
                  export S3_SECRET_KEY="eicS3read"
                '';
              };
          });

    };
}
