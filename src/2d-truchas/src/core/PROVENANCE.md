# 2D simulation provenance

Every JSON-driven 2D simulation should begin its log with compact,
machine-readable prologue records. The records document the executable and
the essential build and run context; they are not a second copy of the input
or a general system inventory.

Each record occupies one line. Its first word identifies the record kind;
the remaining `key=value` fields use double-quoted, backslash-escaped string
values and unquoted numeric values.

The intended fields are:

- caller-supplied code name and Truchas version;
- configured Fortran flags, compiler identity and version, and build host;
- run-start timestamp, run host, and MPI-process count;
- I/O-process CPU model;
- selected library or accelerator versions when they materially affect a
  result.

The original JSON input is copied byte-for-byte to `input.json` in the output
directory. Initialization records will identify that fixed logical name and
its SHA-256 digest rather than embedding the input or serializing its parsed
parameter list.

## Privacy policy

The default prologue must not write the input path, output path, working
directory, build directory, or user name. Those values can expose project or
account information in logs that are later shared. The build and run hosts
are recorded deliberately; users sharing logs should remove them if their
site policy requires it.

The copied input is a separate run artifact, not log content. It may itself
contain paths or other sensitive simulation information, so users running in
a shared output location remain responsible for its access permissions.
