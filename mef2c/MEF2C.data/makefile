build:
	R CMD install --no-test-load ../SingleGeneData
	R CMD install --no-test-load ../MEF2C.data
	R CMD install --no-test-load ../SingleGeneAnalyzer

test:
	R -f MEF2C.data/inst/unitTests/test_MEF2C.data
	R -f SingleGeneAnalyziaer/inst/unitTests/test_SingleGeneAnalyzer.data



