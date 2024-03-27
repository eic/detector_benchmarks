configfile: "config.yaml"

if config["remote"] == "S3":
    from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
    provider = S3RemoteProvider(
        endpoint_url="https://eics3.sdcc.bnl.gov:9000",
        access_key_id=os.environ["S3_ACCESS_KEY"],
        secret_access_key=os.environ["S3_SECRET_KEY"],
    )
    remote_path = lambda path: f"eictest/{path}"
elif config["remote"] == "XRootD":
    from snakemake.remote.XRootD import RemoteProvider as XRootDRemoteProvider
    provider = XRootDRemoteProvider(
        stay_on_remote=False,
    )
    remote_path = lambda path: f"root://dtn-eic.jlab.org//work/eic2/{path}"
else:
    raise ValueError(f"Unexpected config[\"remote\"] = {config['remote']}")

storage = provider.remote # provide Snakemake 8.x interface

include: "Snakefile.common"
