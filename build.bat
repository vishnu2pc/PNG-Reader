cls
@echo off
@if "%1"=="release" (
	set arg1=release
) else (
	set arg1=debug
)

@mkdir bin
@pushd bin
@if %arg1%==release (
	cl ../src/main.cpp /EHsc /O2
	echo Building release
) else (
	cl ../src/main.cpp /EHsc /Zi
	echo Building debug
)
@popd

@pushd src
ctags -R
@popd
