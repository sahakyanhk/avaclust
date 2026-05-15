from setuptools import setup, find_packages

setup(
    name="rnpclust",
    version="0.1.0",
    description="Structural clustering of ribonucleoproteins and RNA-protein interfaces using USalign",
    author="Harutyun Sahakyan",
    license="MIT",
    python_requires=">=3.6",
    install_requires=[
        "numpy",
        "scipy",
        "biopython",
    ],
    packages=find_packages(),
    scripts=[
        "bin/rnpclust",
        "bin/usalign_all_vs_all.bash",
        "bin/align_clusters.bash",
    ],
    entry_points={
        "console_scripts": [
            "rnpclust-hierarchical=bin.hierarchical_cluster:main",
            "rnpclust-setcover=bin.setcover_cluster:main",
            "rnpclust-extract-interface=bin.extract_interface:main",
        ],
    },
)
