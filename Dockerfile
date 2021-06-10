from julia:1.5.1

WORKDIR /code
COPY install.jl .

RUN julia install.jl
