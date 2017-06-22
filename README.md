# linalg
Linear algebra library in Crystal, uses LAPACKE.
- direct access to LAPACK methods
- convenient Matrix(T) class, supports T=Float32, Float64 and Complex.
- high-level interface similar to scipy.linalg or MATLAB. (WIP)

Killing SciPy, one module at a time.

## Installation

1. Install LAPACKE (and dependencies - LAPACK and BLAS). `sudo apt install liblapacke` for Ubuntu, `sudo pacman -S lapacke` for Arch.
2. (for Ubuntu) it seems package doesn't create symlink, so use `sudo ln -s /usr/lib/liblapacke.so.3 /usr/lib/liblapacke.so`

Add this to your application's `shard.yml`:

```yaml
dependencies:
  linalg:
    github: konovod/linalg
```

## Usage

```crystal
require "linalg"
```

TODO: Write examples here.

Check `spec` directory for some examples.

## Development

TODO: Write development instructions here

PRs are welcome)

## Contributing

1. Fork it ( https://github.com/konovod/linalg/fork )
2. Create your feature branch (git checkout -b my-new-feature)
3. Commit your changes (git commit -am 'Add some feature')
4. Push to the branch (git push origin my-new-feature)
5. Create a new Pull Request


## Contributors

- [[konovod]](https://github.com/konovod) Konovod - creator, maintainer
