#!/bin/sh
cd ${0%/*} || exit 1    # Run from this directory

cd ITHACAoutput/Offline2

postProcess -func "minMaxMagnitude(U)"
postProcess -func "patchAverage(name=outlet,U)"
postProcess -func "probes"
postProcess -func "flowRatePatch(name=outlet, phi)"

cd ../Reconstruction

postProcess -func "minMaxMagnitude(U_rec)"
postProcess -func "patchAverage(name=outlet,U_rec)"
postProcess -func "probes"
postProcess -func "flowRatePatch(name=outlet, Phi_rec)"



