# Start with the tidyverse enabled image for R 4.2.0
FROM rocker/tidyverse:4.2.0

# Set environment variables to prevent prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive
# This line will avoid some locale warnings when running as singularity
ENV LC_ALL=C

# Install system dependencies required for R packages
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install R packages: pacman, janitor
RUN R -e "install.packages(c('pacman','janitor'), repos='https://cloud.r-project.org/')"

#specific version of dbplyr
RUN R -e "devtools::install_version('dbplyr', version = '2.3.4', repos = 'https://cloud.r-project.org')"

# Use devtools to install GRETTA from github
RUN R -e "devtools::install_github('MarraLab/GRETTA')"

# Set the default working directory
WORKDIR /usr/local/src/app

### lines below allow people to fire up R (default) or provide an Rscript at runtime
# Default to using R as the entrypoint
ENTRYPOINT ["R"]

# Allow users to pass commands or scripts as arguments
CMD ["--no-save"]

# Set a default command to start R
CMD ["R"]