# MUTLI-STAGE DOCKERFILE for building a Python image
# STAGE 1: builder is the disposable image that makes build products into /opt
# STAGE 2: (default) is the final image that is copied from /opt in the builder image

########################
# STAGE 1: builder     #
########################

# IMPORTANT
# If you change the base image, you will need to update the
# PRE_FETCH_BASE_IMAGE variable in the .gitlab-ci.yml file also.
FROM rocker/r-base:4.2.2 as builder

# Be the root user, explictly
USER root

# Set the top level environment variables
ENV \
    DATA_DIRECTORY="/data" \
    OPT_DIRECTORY="/opt" \
    USER_NAME="admin" \
    USER_DIRECTORY="/home/admin" \
    LC_ALL="en_US.UTF-8" \
    LANG="en_US.UTF-8" 
    
# Set next environment variables that interpolate the top level environment
# variables
ENV \
    USER_BASHRC="${USER_DIRECTORY:?}/.bashrc" \
    USER_BIN_DIRECTORY="${USER_DIRECTORY:?}/.local/bin" \
    SSH_DIR="${USER_DIRECTORY:?}/.ssh" \
    PROJECT_DIRECTORY="${OPT_DIRECTORY:?}/repo" \
    RENV_DIRECTORY="${OPT_DIRECTORY:?}/renv" \
    RENV_PATHS_CACHE="${OPT_DIRECTORY:?}/renv-cache" \
    LOGGING_DIRECTORY="${DATA_DIRECTORY:?}/logs" \
    OTHER_INSTALLS_DIRECTORY="${OPT_DIRECTORY:?}/other-installs" \
    OTHER_INSTALLS_LIB="${OPT_DIRECTORY:?}/other-installs/lib" \
    OTHER_INSTALLS_BIN="${OPT_DIRECTORY:?}/other-installs/bin"

# Set the environment variables for the versions of the software 
ENV \
    PATH="${OTHER_INSTALLS_BIN:?}:$PATH" \
    LD_LIBRARY_PATH="${OTHER_INSTALLS_LIB:?}" \
    RENV_PATHS_LIBRARY="${RENV_DIRECTORY:?}/library"

# Run the commands to:
# - create directories defined in the environment variables
# - create the non-root user
# - put the non-root user in the same groups as the docker user
# - give non-root user ownership of the directories
RUN \
    locale-gen "${LANG:?}" \
    && update-locale LANG="${LANG:?}" \
    && useradd "${USER_NAME}" --shell /bin/bash --create-home --home-dir "${USER_DIRECTORY}" \
    && usermod -a -G docker,staff admin \
    && mkdir -p "${PROJECT_DIRECTORY:?}" "${DATA_DIRECTORY:?}" "${OPT_DIRECTORY:?}" "${RENV_DIRECTORY:?}" "${RENV_PATHS_LIBRARY:?}" "${RENV_PATHS_CACHE:?}" "${OTHER_INSTALLS_LIB:?}" "${OTHER_INSTALLS_BIN:?}" \
    && chown -R "${USER_NAME:?}:${USER_NAME:?}" "${PROJECT_DIRECTORY:?}" "${DATA_DIRECTORY:?}" "${USER_DIRECTORY:?}" "${OPT_DIRECTORY:?}" \
    && chmod -R 755 "${PROJECT_DIRECTORY:?}" "${DATA_DIRECTORY:?}" "${USER_DIRECTORY:?}" "${OPT_DIRECTORY:?}" 

# Install the system dependencies for the R packages but NO R packages
# themselves -- this is done in a later stage.
RUN \
    apt-get update -y \
    && apt-get install -yq --no-install-recommends \
    build-essential \
    apt-transport-https \
    ttf-ubuntu-font-family \
    curl \
    tree \
    cpanminus \
    procps \
    ca-certificates \
    libtasn1-dev \
    nettle-dev \
    libgmp-dev \
    libp11-kit-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libxml2-dev \
    libncurses-dev \
    libudunits2-dev \
    bc \
    && rm -rf /var/lib/apt/lists/*

ENV PATH=/usr/local/bin:$PATH
#######################################
# Install perl v5.38 - required modules
#######################################
# apt-get required libraries like perl, curl and build-essential were added to the code above
# Install cpanminus (cpanm)
RUN curl -L https://cpanmin.us | perl - App::cpanminus

# Install required Perl modules using cpanm
RUN cpanm strict List::MoreUtils List::Util File::Basename Getopt::Long


# Verify installations
RUN perl -Mstrict -e 'print "strict module is available\n"' && \
    perl -MList::MoreUtils -e 'print "List::MoreUtils module is available\n"' && \
    cpanm --version

# Create and CD to the directory where the R packages will be installed
WORKDIR $PROJECT_DIRECTORY

# Copy the R scripts that install the R packages and that parameterize the renv
# The presence of the renv.lock is optional, to achieve idempotency
COPY --chown="${USER_NAME}:${USER_NAME}" [".gitignore", ".Rprofile", "renv.loc[k]", "dependencies.R", "setup_dependencies.R", "./"]
COPY --chown="${USER_NAME}:${USER_NAME}" ["renv/activate.R", "renv/settings.json", "./renv/"]

RUN Rscript setup_dependencies.R --force

# We set the R_LIBS environment variable to the renv library path so that the
# Rscript command will use the renv library path by default (regardless of the
# directory we are in)
ENV R_LIBS="/opt/renv/library/R-4.2/x86_64-pc-linux-gnu"


# Copy the rest of the files into the image (except the .Dockerignore files)
COPY --chown="${USER_NAME}:${USER_NAME}" . .

# # Reapply permissions after all installation and copying is done so the user can
# # manipulate the files if necessary
RUN \
    chown -R "${USER_NAME:?}:${USER_NAME:?}" "${PROJECT_DIRECTORY:?}" "${DATA_DIRECTORY:?}" "${USER_DIRECTORY:?}" "${OPT_DIRECTORY:?}" \
    && chmod -R 755 "${PROJECT_DIRECTORY:?}" "${DATA_DIRECTORY:?}" "${USER_DIRECTORY:?}" "${OPT_DIRECTORY:?}"

# # Switch to the non-root user
USER "${USER_NAME:?}"


# # Build time test to confirm that the package is installed
# # TODO - remove this once CI has a better way to test. We want to check R finds
# # these packages when outside of the repo (the repo is bootstrapped by .Rprofile)
WORKDIR /
RUN R --version && \
    R --slave -e 'packageVersion("dplyr")'

# Explicitly set the working directory (again) and the command to run when the
# container is started
WORKDIR $PROJECT_DIRECTORY
CMD ["/bin/bash"]