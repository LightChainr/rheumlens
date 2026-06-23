from __future__ import annotations

PRIMARY_ISG_15: tuple[str, ...] = (
    "ISG15",
    "IFI6",
    "MX1",
    "OAS1",
    "OAS2",
    "OAS3",
    "IFIT1",
    "IFIT3",
    "IFI44",
    "IFI44L",
    "STAT1",
    "RSAD2",
    "IFITM1",
    "IFITM3",
    "HERC5",
)

AUTHORITATIVE_PRIMARY_RESULTS = {
    "GSE135779": {"scgpt": 0.854, "pca": 0.948, "donor_mean_hvg": 0.898},
    "GSE285773": {"scgpt": 0.950, "pca": 0.988, "donor_mean_hvg": 0.975},
    "GSE174188": {"scgpt": 0.978, "pca": 0.981, "donor_mean_hvg": 0.984},
}

A800_SERVER_PROFILE = {
    "os": "Ubuntu 22.04.4 LTS",
    "gpu": "NVIDIA A800-SXM4-80GB x1",
    "cpu_cores": 144,
    "memory_gib": 1024,
    "torch": "2.5.1+cu124",
    "driver_cuda": "13.0",
    "torch_cuda": "12.4",
    "project_root": "/autodl-fs/data/rheumlens",
    "scratch_root": "/root/autodl-tmp",
}
