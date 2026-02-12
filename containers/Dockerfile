# Using a stable R version (4.3.3 or latest 4.4)
FROM rocker/r-ver:4.3.3

# Install system dependencies for the R packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    libhdf5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    zlib1g-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libglu1-mesa-dev \
    pandoc \
    build-essential \
    gfortran \
    libblas-dev \
    liblapack-dev \
    cmake \
    texlive-latex-base \
    texlive-latex-recommended \
    texlive-fonts-recommended \
    texlive-latex-extra \
    texlive-fonts-extra \
    && rm -rf /var/lib/apt/lists/*


# Install BiocManager for DESeq2
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org/')"
RUN R -e "BiocManager::install(c('DESeq2','rhdf5'))"

# Install CRAN packages
RUN R -e "install.packages(c(\
    'future', 'Matrix', 'Seurat', 'dplyr', 'tidyr', \
    'patchwork', 'ggplot2', 'rmarkdown', \
    'lme4', 'emmeans', 'argparser','ggpubr', \
    'jsonlite', 'RColorBrewer', 'ggtext', 'tinytex'\
    ), dependencies=TRUE, repos='https://cloud.r-project.org/')"

# Initialize tinytex
RUN R -e "tinytex::install_tinytex(force = TRUE)"

CMD ["R"]