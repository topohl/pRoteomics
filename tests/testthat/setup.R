# No test may append to the reviewer ledger.
#
# results/reviewer_audit/input_resolution_audit.csv is production evidence:
# reviewers read it to see how every scientific input was resolved. In Phase
# 6I.2 a unit test put 80 rows into it naming a drive letter that does not
# exist, and they had to be found and removed again. The general problem is
# that a test and the pipeline resolve inputs through the same function, and
# that function records provenance to one fixed path.
#
# Routing the destination fixes it for the whole suite in one place. Doing it
# here rather than per test file matters for two reasons: a file that forgets
# is silently back to writing production evidence, and testthat sources this
# before every test file, including when a single file is run on its own.
#
# The override is an environment variable on purpose. Eleven test files spawn
# real analysis scripts with system2(), and a child process inherits the
# environment - an in-process fixture would not reach them.
#
# Isolation is by routing, not by a dry-run guard in the writer. The appender
# has to stay unconditional: test-preprocessing-writer-namespace.R asserts it
# contains no dry-run guard, because that is exactly why
# build_module_score_metadata is classified PATH_VERIFIED_STRUCTURALLY_ONLY.
#
# A test that wants to assert on what was written should call
# local_input_resolution_audit() to point at its own file; that override takes
# precedence over this one and is restored on exit.

local({
  ledger <- file.path(tempdir(), "testthat-input-resolution-audit.csv")
  Sys.setenv(PROTEOMICS_INPUT_RESOLUTION_AUDIT = ledger)
  withr::defer(Sys.unsetenv("PROTEOMICS_INPUT_RESOLUTION_AUDIT"),
               testthat::teardown_env())
})
