# Installing ZarzaETC

At the present time, the only tested and officially supported platforms are macOS and GNU/Linux. There are no prebuilt binaries, so you will need to build ZarzaETC from source.

## Requirements

Before building ZarzaETC, you must ensure that you have the following installed:

* [Rust](https://rust-lang.org/es/tools/install/) with Cargo, version 1.70 or newer
* [Git](https://www.git-scm.com/) or later

To install Rust, you may use `rustup`, which also installs `cargo`.
On Linux and macOS systems, this is done as follows:

```bash
$ curl https://sh.rustup.rs -sSf | sh
```

## Download the code

You can download the code of the TARSIS ETC with:

```bash
$ git clone https://github.com/HARMONI-CAB/ZarzaETC
```

## Install and compile

ZarzaETC is a typical rust based program for which `cargo` is necessary to compile. Just run, in order:

```bash
$ cd ZarzaETC
$ cargo build
```

If the installation and compilation was succesfull, you can now run ZarzaETC CLI as:

```bash
$ zarza_etc -h
```

to display usage information.
