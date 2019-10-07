# specfemx_benchmark_assist

Python program to calculate magnitude of displacement results from SPECFEMX. It integrates over the domain using gll quadrature.
It can be used to compare two results (for benchmarking) or calculate the net displacement in a domain.
When comparing results one must be resampled so they are on exactly the same integration points.

An example is included to demonstrate how the package should be used -- a comparison of SPECFEMX and Okada solutions. A file should be created that follows the format of the compare_pkada_hn_cmt.py file. The following mesh information is needed (see the Okada folder): mesh coordinates and connectivity, and a list of elements that comprise the fault or contain CMT sources (these elements are not included in the integration). The two solutions to compare should be provided in .csv format.

The difference between the two solutions is calculated by integrating over the difference at the GLL points using the spectral element quadrature and dividing the result by the integrated deformation from the reference solution (generally Okada or Pylith).
