# syntax=docker/dockerfile:1

# Debian release shared by the builder and runtime stages. Keep them in sync.
ARG DEBIAN_RELEASE=trixie

# ---- Builder: compile the apoe-typer CLI --------------------------------
FROM rust:1-slim-${DEBIAN_RELEASE} AS builder
WORKDIR /src
COPY . .
RUN cargo build --release --locked --bin apoe-typer

# ---- Test: run the unit and integration tests (not the default target) --
FROM builder AS test
RUN cargo test --release --locked

# ---- Runtime: minimal image with only the binary ------------------------
FROM debian:${DEBIAN_RELEASE}-slim AS runtime

LABEL org.opencontainers.image.title="apoe-typer" \
      org.opencontainers.image.description="APOE genotyper (rs429358 / rs7412) for VCF and gVCF files" \
      org.opencontainers.image.version="0.1.0"

COPY --from=builder /src/target/release/apoe-typer /usr/local/bin/apoe-typer

# Unprivileged user with a fixed numeric UID so bind-mounted files can be matched.
RUN useradd --uid 10001 --user-group --no-create-home --shell /usr/sbin/nologin apoe \
    && mkdir -p /data \
    && chown apoe:apoe /data

USER apoe:apoe
WORKDIR /data
ENTRYPOINT ["apoe-typer"]
