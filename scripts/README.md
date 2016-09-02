Scripts
------------------------------------------------------------------------------
This folder contains various scripts that may be useful for users.

### Build Scripts
#### [build](./build)
Builds both Truchas and the libraries that it depends on.  Usage:

    build <install_location> <compiler> <mode> <build>

* `<install_location>`: Directory to install Truchas and required libaries
* `<compiler>`: Compiler to use. Options: `intel`, `nag`, `gcc`
* `<mode>`: Mode in which Truchas runs. Options: `serial`, `parallel`
* `<build>`: Build type. Options: `dbg` for debug, `opt` for release, and
  `check` for check
* See [build_tpl](./build_tpl) and [build_truchas](./build_truchas) for more
  details.

#### [build_tpl](./build_tpl)
Builds only the libraries required for Truchas.  Usage:

    build_tpl <install_location> <compiler> <mode>

* `<install_location>`, `<compiler>`, `<mode>` are the same as above.
* See [tpl/config](../tpl/config) for more details

#### [build_truchas](./build_truchas)
Builds only Truchas; required libraries must exist prior to calling this
script.  Usage:

    build_truchas <install_location> <tpl_location> <compiler> <mode> <build>

* `<install_location>`, `<compiler>`, `<mode>`, `<build>` are the same as
    above.
* See [config](../config) for more details.

#### [test_truchas](./test_truchas)
Runs regression tests for Truchas.  Usage:

    test_truchas <truchas_location> <num_of_proc>

* `<truchas_location>`: Location where Truchas was installed with
    [build](./build) or [build_truchas](./build_truchas).
* `<num_of_proc>`: Number of processes to use for testing
