FROM ghcr.io/prefix-dev/pixi:0.82.0 AS build

WORKDIR /app
COPY pixi.lock pyproject.toml LICENSE /app/
COPY src /app/src
RUN pixi install --locked \
    && pixi shell-hook -e default -s bash > /app/activate.sh \
    && /app/.pixi/envs/default/bin/python -m pip uninstall -y preface \
    && /app/.pixi/envs/default/bin/python -m pip install --no-deps .

FROM ubuntu:26.04 AS production

WORKDIR /app
COPY --from=build /app/.pixi/envs/default /app/.pixi/envs/default
COPY --from=build /app/activate.sh /app/activate.sh

RUN printf '#!/bin/bash\n. /app/activate.sh\nexec /app/.pixi/envs/default/bin/preface "$@"\n' \
    > /usr/local/bin/preface \
    && chmod 0755 /usr/local/bin/preface

CMD ["bash"]
