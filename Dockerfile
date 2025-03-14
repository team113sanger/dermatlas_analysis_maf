# MUTLI-STAGE DOCKERFILE for building a R image
# STAGE 1: builder is the disposable image that makes build products into /opt

########################
# STAGE 1: builder     #
########################

# IMPORTANT
# If you change the base image, you will need to update the
# PRE_FETCH_BASE_IMAGE variable in the .gitlab-ci.yml file also.
FROM rocker/tidyverse:4.2.2 as builder

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
    LOGGING_DIRECTORY="${DATA_DIRECTORY:?}/logs" 
# Set the environment variables for the versions of the software 
ENV \
    RENV_PATHS_LIBRARY="${RENV_DIRECTORY:?}/library"

# Run the commands to:
# - create directories defined in the environment variables
# - create the non-root user
# - put the non-root user in the same groups as the docker user
# - give non-root user ownership of the directories

# Be the root user, explictly
USER root

RUN \
    locale-gen "${LANG:?}" \
    && update-locale LANG="${LANG:?}" \
    && useradd "${USER_NAME}" --shell /bin/bash --create-home --home-dir "${USER_DIRECTORY}" \
    && groupadd -r docker \
    && usermod -a -G docker,staff admin \
    && mkdir -p "${PROJECT_DIRECTORY:?}" "${DATA_DIRECTORY:?}" "${OPT_DIRECTORY:?}" "${RENV_DIRECTORY:?}" "${RENV_PATHS_LIBRARY:?}" "${RENV_PATHS_CACHE:?}" \
    && chown -R "${USER_NAME:?}:${USER_NAME:?}" "${PROJECT_DIRECTORY:?}" "${DATA_DIRECTORY:?}" "${USER_DIRECTORY:?}" "${OPT_DIRECTORY:?}" \
    && chmod -R 755 "${PROJECT_DIRECTORY:?}" "${DATA_DIRECTORY:?}" "${USER_DIRECTORY:?}" "${OPT_DIRECTORY:?}" 

# Install the system dependencies for the R packages but NO R packages
# themselves -- this is done in a later stage.
RUN \
    apt-get update -y \
    && apt-get install -yq --no-install-recommends \
    build-essential \
    apt-transport-https \
    curl \
    tree \
    cpanminus \
    procps \
    ca-certificates \
    libtasn1-dev \
    nettle-dev \
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
COPY --chown="${USER_NAME}:${USER_NAME}" [".gitignore", "renv.loc[k]", "./"]
# COPY --chown="${USER_NAME}:${USER_NAME}" ["renv/activate.R", "renv/settings.json", "./renv/"]

ENV SYSTEM_LIBS="/usr/local/lib/R/site-library"
RUN R --slave -e "install.packages('renv')"
RUN R --slave -e "renv::init(settings = list(external.libraries = '${SYSTEM_LIBS}'))"

# We set the R_LIBS environment variable to the renv library path so that the
# Rscript command will use the renv library path by default (regardless of the
# directory we are in)
RUN R -e "renv::restore(prompt= FALSE)"

# RUN Rscript setup_dependencies.R --force
ENV R_LIBS="/opt/renv/library/R-4.2/x86_64-pc-linux-gnu:${SYSTEM_LIBS:?}:/opt/renv-cache/v5/R-4.2/x86_64-pc-linux-gnu"



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
    R --slave -e "packageVersion('vroom')" && \
    R --slave -e "packageVersion('cowplot')"

    # Explicitly set the working directory (again) and the command to run when the
# container is started
WORKDIR $PROJECT_DIRECTORY
CMD ["/bin/bash"]