#!/bin/sh
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
state="$here/seven_line_stage33_raw.json"

jq -r '.charts[] | select(.id==212) |
  "ring r=0,(x0,x2,x3),dp; poly f=" + .polynomial + ";
   ideal C1=x2*x3^2+x3^3-4*x2*x3-4*x3^2+4*x2+3*x3,
     2*x0*x3^2+x2*x3^2+x3^3-2*x0*x3-4*x2*x3-4*x3^2+4*x2+x3+4,
     2*x0*x2-x2*x3-x3^2+2*x2+2*x3+3;
   ideal G1=std(C1);
   ideal R1=reduce(f,G1),reduce(diff(f,x0),G1),reduce(diff(f,x2),G1),reduce(diff(f,x3),G1);
   print(\"CHART_212_C1_DIM=\"+string(dim(G1))); print(R1);
   ideal C2=x2+x3^2,x0*x3^2-x3-1;
   ideal G2=std(C2);
   ideal R2=reduce(f,G2),reduce(diff(f,x0),G2),reduce(diff(f,x2),G2),reduce(diff(f,x3),G2);
   print(\"CHART_212_C2_DIM=\"+string(dim(G2))); print(R2); quit;"' "$state" | Singular -q

jq -r '.charts[] | select(.id==213) |
  "ring r=0,(x0,x2,x3),dp; poly f=" + .polynomial + ";
   ideal C1=x2*x3-2*x2+x3-3,
     x0*x3^3-4*x0*x3^2+3*x0*x3+x3^2-4*x3+4,
     -x0*x3^2+2*x0*x2+2*x0*x3+3*x0-x3+2;
   ideal G1=std(C1);
   ideal R1=reduce(f,G1),reduce(diff(f,x0),G1),reduce(diff(f,x2),G1),reduce(diff(f,x3),G1);
   print(\"CHART_213_C1_DIM=\"+string(dim(G1))); print(R1);
   ideal C2=x0*x3^2+1,x2+x3+1;
   ideal G2=std(C2);
   ideal R2=reduce(f,G2),reduce(diff(f,x0),G2),reduce(diff(f,x2),G2),reduce(diff(f,x3),G2);
   print(\"CHART_213_C2_DIM=\"+string(dim(G2))); print(R2); quit;"' "$state" | Singular -q
