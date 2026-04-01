import { useEffect, useRef, useState } from 'react'
import { Card, Button } from '@/components/ui'

export interface ViewerProps {
  trajectoryPath?: string
  topologyPath?: string
  initialFrame?: number
  rmsdDataPath?: string
  bindingSite?: { center: [number, number, number]; radius: number }
}

interface TrajectoryPlayerState {
  isPlaying: boolean
  currentFrame: number
  totalFrames: number
  playbackSpeed: number
  rmsdHistory: number[]
  time_ps: number
}

interface BindingSite {
  center: [number, number, number]
  radius: number
}

const VIEW_STYLES = [
  { id: 'stick', label: 'Stick' },
  { id: 'ball', label: 'Ball & Stick' },
  { id: 'cartoon', label: 'Cartoon' },
  { id: 'surface', label: 'Surface' },
  { id: 'sphere', label: 'Spacefill' },
]

const COLOR_SCHEMES = [
  { id: 'element', label: 'Element' },
  { id: 'sstruc', label: 'Secondary Structure' },
  { id: 'chain', label: 'Chain' },
  { id: 'bfactor', label: 'B-Factor' },
  { id: 'chainHetatm', label: 'Chains (Colorful)' },
]

// Sample protein structure with multiple chains for colorful display
const SAMPLE_PROTEIN = `
HEADER    IMMUNE SYSTEM                           15-APR-98   1HIA
TITLE     CRYSTAL STRUCTURE OF HUMAN HEAT-LABILE ENTEROTOXIN
EXPDTA    X-RAY DIFFRACTION
AUTHOR    M.E. NEMETH,J.M.BADANO,A.C.WALLIS,B.LOCHNER,A.HOLZENBURG
HETATM    1  N   MET A   1       9.287  11.614  -1.060  1.00 28.12           N
HETATM    2  CA  MET A   1      10.286  11.038  -1.989  1.00 27.95           C
HETATM    3  C   MET A   1      11.669  11.451  -1.661  1.00 28.26           C
HETATM    4  O   MET A   1      12.039  12.619  -1.589  1.00 30.27           O
HETATM    5  N   ILE A   2      12.516  10.421  -1.472  1.00 27.21           N
HETATM    6  CA  ILE A   2      13.922  10.668  -1.148  1.00 27.26           C
HETATM    7  C   ILE A   2      14.802  10.296  -2.322  1.00 27.85           C
HETATM    8  O   ILE A   2      15.995  10.583  -2.399  1.00 29.59           O
HETATM    9  N   GLN A   3      14.261   9.687  -3.311  1.00 26.84           N
HETATM   10  CA  GLN A   3      14.996   9.257  -4.500  1.00 27.01           C
HETATM   11  C   GLN A   3      14.103   8.443  -5.420  1.00 28.07           C
HETATM   12  O   GLN A   3      13.020   7.954  -5.086  1.00 28.77           O
HETATM   13  N   VAL A   4      14.578   8.279  -6.640  1.00 26.58           N
HETATM   14  CA  VAL A   4      13.816   7.510  -7.636  1.00 26.74           C
HETATM   15  C   VAL A   4      12.339   7.896  -7.742  1.00 27.03           C
HETATM   16  O   VAL A   4      11.513   7.663  -6.843  1.00 27.72           O
HETATM   17  N   GLU A   5      12.000   8.508  -8.879  1.00 26.07           N
HETATM   18  CA  GLU A   5      10.630   8.932  -9.112  1.00 26.23           C
HETATM   19  C   GLU A   5      10.523  10.453  -9.246  1.00 27.04           C
HETATM   20  O   GLU A   5       9.424  11.008  -9.378  1.00 28.03           O
HETATM   21  N   TYR A   6       9.876  11.657 -10.220  1.00 25.97           N
HETATM   22  CA  TYR A   6       9.847  13.112 -10.446  1.00 26.46           C
HETATM   23  C   TYR A   6       8.484  13.715 -10.165  1.00 27.33           C
HETATM   24  O   TYR A   6       7.453  13.062  -9.997  1.00 27.68           O
HETATM   25  N   CYS A   7       8.484  15.049 -10.107  1.00 26.46           N
HETATM   26  CA  CYS A   7       7.236  15.757 -9.845  1.00 27.05           C
HETATM   27  C   CYS A   7       7.390  17.269 -9.673  1.00 27.96           C
HETATM   28  O   CYS A   7       8.511  17.810  -9.730  1.00 29.03           O
HETATM   29  N   GLY A   8       6.251  17.960  -9.483  1.00 26.73           N
HETATM   30  CA  GLY A   8       6.291  19.410  -9.296  1.00 26.85           C
HETATM   31  C   GLY A   8       7.007  20.124 -10.412  1.00 27.35           C
HETATM   32  O   GLY A   8       7.717  19.489 -11.226  1.00 27.81           O
HETATM   33  N   GLU A   9       6.847  21.437 -10.500  1.00 26.51           N
HETATM   34  CA  GLU A   9       7.486  22.239 -11.545  1.00 27.08           C
HETATM   35  C   GLU A   9       8.986  22.001 -11.640  1.00 27.91           C
HETATM   36  O   GLU A   9       9.585  21.047 -11.107  1.00 28.92           O
HETATM   37  N   THR A  10       9.634  22.902 -12.310  1.00 27.29           N
HETATM   38  CA  THR A  10      11.084  22.769 -12.452  1.00 27.72           C
HETATM   39  C   THR A  10      11.842  23.057 -11.165  1.00 28.32           C
HETATM   40  O   THR A  10      11.225  23.320 -10.120  1.00 29.04           O
HETATM   41  N   GLY A  11      13.183  23.026 -11.251  1.00 27.44           N
HETATM   42  CA  GLY A  11      14.013  23.310 -10.083  1.00 27.68           C
HETATM   43  C   GLY A  11      13.689  22.362  -8.939  1.00 28.04           C
HETATM   44  O   GLY A  11      14.587  22.020  -8.152  1.00 28.73           O
HETATM   45  N   ASN A  12      12.412  21.934  -8.847  1.00 27.27           N
HETATM   46  CA  ASN A  12      11.975  21.019  -7.794  1.00 27.42           C
HETATM   47  C   ASN A  12      12.231  21.578  -6.409  1.00 27.82           C
HETATM   48  O   ASN A  12      12.717  22.707  -6.241  1.00 28.53           O
HETATM   49  N   ALA A  13      11.932  20.780  -5.372  1.00 27.07           N
HETATM   50  CA  ALA A  13      12.133  21.230  -4.002  1.00 27.03           C
HETATM   51  C   ALA A  13      10.828  21.234  -3.219  1.00 27.27           C
HETATM   52  O   ALA A  13      10.853  21.522  -2.011  1.00 28.05           O
HETATM   53  N   ALA A  14       9.695  20.956  -3.872  1.00 26.52           N
HETATM   54  CA  ALA A  14       8.434  20.910  -3.181  1.00 26.45           C
HETATM   55  C   ALA A  14       7.290  20.857  -4.182  1.00 26.69           C
HETATM   56  O   ALA A  14       7.490  20.706  -5.398  1.00 27.37           O
HETATM   57  N   VAL A  15       6.062  20.957  -3.667  1.00 25.88           N
HETATM   58  CA  VAL A  15       4.869  20.917  -4.505  1.00 26.16           C
HETATM   59  C   VAL A  15       4.966  21.891  -5.673  1.00 26.87           C
HETATM   60  O   VAL A  15       5.996  22.579  -5.831  1.00 27.85           O
HETATM   61  N   ASP A  16       3.910  21.950  -6.500  1.00 25.87           N
HETATM   62  CA  ASP A  16       3.913  22.857  -7.653  1.00 26.33           C
HETATM   63  C   ASP A  16       4.350  24.286  -7.313  1.00 27.19           C
HETATM   64  O   ASP A  16       4.775  24.587  -6.182  1.00 28.05           O
HETATM   65  N   ALA A  17       4.253  25.188  -8.327  1.00 26.59           N
HETATM   66  CA  ALA A  17       4.647  26.574  -8.137  1.00 26.95           C
HETATM   67  C   ALA A  17       5.892  26.672  -7.282  1.00 27.60           C
HETATM   68  O   ALA A  17       6.275  27.746  -6.807  1.00 28.38           O
HETATM   69  N   SER A  18       6.528  25.542  -7.065  1.00 26.85           N
HETATM   70  CA  SER A  18       7.733  25.496  -6.260  1.00 27.24           C
HETATM   71  C   SER A  18       7.477  25.202  -4.789  1.00 27.71           C
HETATM   72  O   SER A  18       6.338  24.888  -4.407  1.00 28.38           O
HETATM   73  N   LEU A  19       8.517  25.297  -3.941  1.00 26.96           N
HETATM   74  CA  LEU A  19       8.364  25.038  -2.513  1.00 27.18           C
HETATM   75  C   LEU A  19       9.700  24.841  -1.825  1.00 27.70           C
HETATM   76  O   LEU A  19      10.748  24.882  -2.467  1.00 28.39           O
HETATM   77  N   GLN A  20       9.659  24.650  -0.516  1.00 26.99           N
HETATM   78  CA  GLN A  20      10.896  24.452   0.228  1.00 27.35           C
HETATM   79  C   GLN A  20      11.644  25.754   0.495  1.00 28.09           C
HETATM   80  O   GLN A  20      12.850  25.863   0.264  1.00 28.91           O
HETATM   81  N   THR A  21      10.951  26.807   0.964  1.00 27.35           N
HETATM   82  CA  THR A  21      11.600  28.097   1.241  1.00 27.63           C
HETATM   83  C   THR A  21      11.847  28.936   0.001  1.00 28.19           C
HETATM   84  O   THR A  21      10.927  29.236  -0.764  1.00 28.78           O
HETATM   85  N   SER A  22      13.104  29.326  -0.255  1.00 27.62           N
HETATM   86  CA  SER A  22      13.475  30.135  -1.422  1.00 27.88           C
HETATM   87  C   SER A  22      12.310  30.915  -2.017  1.00 28.24           C
HETATM   88  O   SER A  22      12.513  31.831  -2.832  1.00 29.00           O
HETATM   89  N   THR A  23      11.089  30.585  -1.627  1.00 27.55           N
HETATM   90  CA  THR A  23       9.893  31.255  -2.108  1.00 27.72           C
HETATM   91  C   THR A  23       9.975  32.754  -1.861  1.00 28.25           C
HETATM   92  O   THR A  23       9.020  33.489  -2.138  1.00 28.93           O
HETATM   93  N   ARG A  24      11.142  33.246  -1.344  1.00 27.68           N
HETATM   94  CA  ARG A  24      11.330  34.681  -1.065  1.00 27.99           C
HETATM   95  C   ARG A  24      11.201  34.953   0.427  1.00 28.47           C
HETATM   96  O   ARG A  24      12.162  34.774   1.181  1.00 29.16           O
HETATM   97  N   THR B  25      14.287  11.614  -1.060  1.00 28.12           N
HETATM   98  CA  THR B  25      15.286  11.038  -1.989  1.00 27.95           C
HETATM   99  C   THR B  25      16.669  11.451  -1.661  1.00 28.26           C
HETATM  100  O   THR B  25      17.039  12.619  -1.589  1.00 30.27           O
HETATM  101  N   ILE B  26      17.516  10.421  -1.472  1.00 27.21           N
HETATM  102  CA  ILE B  26      18.922  10.668  -1.148  1.00 27.26           C
HETATM  103  C   ILE B  26      19.802  10.296  -2.322  1.00 27.85           C
HETATM  104  O   ILE B  26      20.995  10.583  -2.399  1.00 29.59           O
HETATM  105  N   GLN B  27      19.261   9.687  -3.311  1.00 26.84           N
HETATM  106  CA  GLN B  27      19.996   9.257  -4.500  1.00 27.01           C
HETATM  107  C   GLN B  27      19.103   8.443  -5.420  1.00 28.07           C
HETATM  108  O   GLN B  27      18.020   7.954  -5.086  1.00 28.77           O
HETATM  109  N   VAL B  28      19.578   8.279  -6.640  1.00 26.58           N
HETATM  110  CA  VAL B  28      18.816   7.510  -7.636  1.00 26.74           C
HETATM  111  C   VAL B  28      17.339   7.896  -7.742  1.00 27.03           C
HETATM  112  O   VAL B  28      16.513   7.663  -6.843  1.00 27.72           O
HETATM  113  N   GLU B  29      17.000   8.508  -8.879  1.00 26.07           N
HETATM  114  CA  GLU B  29      15.630   8.932  -9.112  1.00 26.23           C
HETATM  115  C   GLU B  29      15.523  10.453  -9.246  1.00 27.04           C
HETATM  116  O   GLU B  29      14.424  11.008  -9.378  1.00 28.03           O
HETATM  117  N   TYR B  30      14.876  11.657 -10.220  1.00 25.97           N
HETATM  118  CA  TYR B  30      14.847  13.112 -10.446  1.00 26.46           C
HETATM  119  C   TYR B  30      13.484  13.715 -10.165  1.00 27.33           C
HETATM  120  O   TYR B  30      12.453  13.062  -9.997  1.00 27.68           O
HETATM  121  N   CYS B  31      13.484  15.049 -10.107  1.00 26.46           N
HETATM  122  CA  CYS B  31      12.236  15.757  -9.845  1.00 27.05           C
HETATM  123  C   CYS B  31      12.390  17.269  -9.673  1.00 27.96           C
HETATM  124  O   CYS B  31      13.511  17.810  -9.730  1.00 29.03           O
HETATM  125  N   GLY B  32      11.251  17.960  -9.483  1.00 26.73           N
HETATM  126  CA  GLY B  32      11.291  19.410  -9.296  1.00 26.85           C
HETATM  127  C   GLY B  32      12.007  20.124 -10.412  1.00 27.35           C
HETATM  128  O   GLY B  32      12.717  19.489 -11.226  1.00 27.81           O
HETATM  129  N   GLU B  33      11.847  21.437 -10.500  1.00 26.51           N
HETATM  130  CA  GLU B  33      12.486  22.239 -11.545  1.00 27.08           C
HETATM  131  C   GLU B  33      13.986  22.001 -11.640  1.00 27.91           C
HETATM  132  O   GLU B  33      14.585  21.047 -11.107  1.00 28.92           O
HETATM  133  N   THR B  34      14.634  22.902 -12.310  1.00 27.29           N
HETATM  134  CA  THR B  34      16.084  22.769 -12.452  1.00 27.72           C
HETATM  135  C   THR B  34      16.842  23.057 -11.165  1.00 28.32           C
HETATM  136  O   THR B  34      16.225  23.320 -10.120  1.00 29.04           O
HETATM  137  N   GLY B  35      18.183  23.026 -11.251  1.00 27.44           N
HETATM  138  CA  GLY B  35      19.013  23.310 -10.083  1.00 27.68           C
HETATM  139  C   GLY B  35      18.689  22.362  -8.939  1.00 28.04           C
HETATM  140  O   GLY B  35      19.587  22.020  -8.152  1.00 28.73           O
HETATM  141  N   ASN B  36      17.412  21.934  -8.847  1.00 27.27           N
HETATM  142  CA  ASN B  36      16.975  21.019  -7.794  1.00 27.42           C
HETATM  143  C   ASN B  36      17.231  21.578  -6.409  1.00 27.82           C
HETATM  144  O   ASN B  36      17.717  22.707  -6.241  1.00 28.53           O
HETATM  145  N   ALA B  37      16.932  20.780  -5.372  1.00 27.07           N
HETATM  146  CA  ALA B  37      17.133  21.230  -4.002  1.00 27.03           C
HETATM  147  C   ALA B  37      15.828  21.234  -3.219  1.00 27.27           C
HETATM  148  O   ALA B  37      15.853  21.522  -2.011  1.00 28.05           O
HETATM  149  N   ALA B  38      14.695  20.956  -3.872  1.00 26.52           N
HETATM  150  CA  ALA B  38      13.434  20.910  -3.181  1.00 26.45           C
HETATM  151  C   ALA B  38      12.290  20.857  -4.182  1.00 26.69           C
HETATM  152  O   ALA B  38      12.490  20.706  -5.398  1.00 27.37           O
HETATM  153  N   VAL B  39      11.062  20.957  -3.667  1.00 25.88           N
HETATM  154  CA  VAL B  39       9.869  20.917  -4.505  1.00 26.16           C
HETATM  155  C   VAL B  39       9.966  21.891  -5.673  1.00 26.87           C
HETATM  156  O   VAL B  39      10.996  22.579  -5.831  1.00 27.85           O
HETATM  157  N   ASP B  40       8.910  21.950  -6.500  1.00 25.87           N
HETATM  158  CA  ASP B  40       8.913  22.857  -7.653  1.00 26.33           C
HETATM  159  C   ASP B  40       9.350  24.286  -7.313  1.00 27.19           C
HETATM  160  O   ASP B  40       9.775  24.587  -6.182  1.00 28.05           O
HETATM  161  N   ALA B  41       9.253  25.188  -8.327  1.00 26.59           N
HETATM  162  CA  ALA B  41       9.647  26.574  -8.137  1.00 26.95           C
HETATM  163  C   ALA B  41      10.892  26.672  -7.282  1.00 27.60           C
HETATM  164  O   ALA B  41      11.275  27.746  -6.807  1.00 28.38           O
HETATM  165  N   SER B  42      11.528  25.542  -7.065  1.00 26.85           N
HETATM  166  CA  SER B  42      12.733  25.496  -6.260  1.00 27.24           C
HETATM  167  C   SER B  42      12.477  25.202  -4.789  1.00 27.71           C
HETATM  168  O   SER B  42      11.338  24.888  -4.407  1.00 28.38           O
HETATM  169  N   LEU B  43      13.517  25.297  -3.941  1.00 26.96           N
HETATM  170  CA  LEU B  43      13.364  25.038  -2.513  1.00 27.18           C
HETATM  171  C   LEU B  43      14.700  24.841  -1.825  1.00 27.70           C
HETATM  172  O   LEU B  43      15.748  24.882  -2.467  1.00 28.39           O
HETATM  173  N   GLN B  44      14.659  24.650  -0.516  1.00 26.99           N
HETATM  174  CA  GLN B  44      15.896  24.452   0.228  1.00 27.35           C
HETATM  175  C   GLN B  44      16.644  25.754   0.495  1.00 28.09           C
HETATM  176  O   GLN B  44      17.850  25.863   0.264  1.00 28.91           O
HETATM  177  N   THR B  45      15.951  26.807   0.964  1.00 27.35           N
HETATM  178  CA  THR B  45      16.600  28.097   1.241  1.00 27.63           C
HETATM  179  C   THR B  45      16.847  28.936   0.001  1.00 28.19           C
HETATM  180  O   THR B  45      15.927  29.236  -0.764  1.00 28.78           O
HETATM  181  N   SER B  46      18.104  29.326  -0.255  1.00 27.62           N
HETATM  182  CA  SER B  46      18.475  30.135  -1.422  1.00 27.88           C
HETATM  183  C   SER B  46      17.310  30.915  -2.017  1.00 28.24           C
HETATM  184  O   SER B  46      17.513  31.831  -2.832  1.00 29.00           O
HETATM  185  N   THR B  47      16.089  30.585  -1.627  1.00 27.55           N
HETATM  186  CA  THR B  47      14.893  31.255  -2.108  1.00 27.72           C
HETATM  187  C   THR B  47      14.975  32.754  -1.861  1.00 28.25           C
HETATM  188  O   THR B  47      14.020  33.489  -2.138  1.00 28.93           O
HETATM  189  N   ARG B  48      16.142  33.246  -1.344  1.00 27.68           N
HETATM  190  CA  ARG B  48      16.330  34.681  -1.065  1.00 27.99           C
HETATM  191  C   ARG B  48      16.201  34.953   0.427  1.00 28.47           C
HETATM  192  O   ARG B  48      17.162  34.774   1.181  1.00 29.16           O
CONECT    1    2
CONECT    2    1    3
CONECT    3    2    4    5
CONECT    4    3
CONECT    5    3    6
CONECT    6    5    7
CONECT    7    6    8    9
CONECT    8    7
CONECT    9    7   10
CONECT   10    9   11
CONECT   11   10   12   13
CONECT   12   11
CONECT   13   11   14
CONECT   14   13   15
CONECT   15   14   16   17
CONECT   16   15
CONECT   17   15   18
CONECT   18   17   19
CONECT   19   18   20   21
CONECT   20   19
CONECT   21   19   22
CONECT   22   21   23
CONECT   23   22   24   25
CONECT   24   23
CONECT   25   23   26
CONECT   26   25   27
CONECT   27   26   28   29
CONECT   28   27
CONECT   29   27   30
CONECT   30   29   31
CONECT   31   30   32   33
CONECT   32   31
CONECT   33   31   34
CONECT   34   33   35
CONECT   35   34   36   37
CONECT   36   35
CONECT   37   35   38
CONECT   38   37   39
CONECT   39   38   40   41
CONECT   40   39
CONECT   41   39   42
CONECT   42   41   43
CONECT   43   42   44   45
CONECT   44   43
CONECT   45   43   46
CONECT   46   45   47
CONECT   47   46   48   49
CONECT   48   47
CONECT   49   47   50
CONECT   50   49   51
CONECT   51   50   52   53
CONECT   52   51
CONECT   53   51   54
CONECT   54   53   55
CONECT   55   54   56   57
CONECT   56   55
CONECT   57   55   58
CONECT   58   57   59
CONECT   59   58   60   61
CONECT   60   59
CONECT   61   59   62
CONECT   62   61   63
CONECT   63   62   64   65
CONECT   64   63
CONECT   65   63   66
CONECT   66   65   67
CONECT   67   66   68   69
CONECT   68   67
CONECT   69   67   70
CONECT   70   69   71
CONECT   71   70   72   73
CONECT   72   71
CONECT   73   71   74
CONECT   74   73   75
CONECT   75   74   76   77
CONECT   76   75
CONECT   77   75   78
CONECT   78   77   79
CONECT   79   78   80   81
CONECT   80   79
CONECT   81   79   82
CONECT   82   81   83
CONECT   83   82   84   85
CONECT   84   83
CONECT   85   83   86
CONECT   86   85   87
CONECT   87   86   88   89
CONECT   88   87
CONECT   89   87   90
CONECT   90   89   91
CONECT   91   90   92   93
CONECT   92   91
CONECT   93   91   94
CONECT   94   93   95
CONECT   95   94   96
CONECT   96   95
END
`

export function Viewer() {
  const containerRef = useRef<HTMLDivElement>(null)
  const viewerRef = useRef<any>(null)
  const [style, setStyle] = useState('cartoon')
  const [colorScheme, setColorScheme] = useState('chain')
  const [loaded, setLoaded] = useState(false)
  const [loadedFileName, setLoadedFileName] = useState<string | null>(null)
  const [_trajectoryPath, setTrajectoryPath] = useState<string | null>(null)
  const [_topologyPath, setTopologyPath] = useState<string | null>(null)
  const [bindingSite, _setBindingSite] = useState<BindingSite | null>(null)
  const [trajState, setTrajState] = useState<TrajectoryPlayerState>({
    isPlaying: false,
    currentFrame: 0,
    totalFrames: 0,
    playbackSpeed: 15,
    rmsdHistory: [],
    time_ps: 0,
  })
  const [showTrajectoryControls, setShowTrajectoryControls] = useState(false)
  const [showSurface, setShowSurface] = useState(false)
  const [_activePanel, _setActivePanel] = useState<'style' | 'trajectory' | 'analysis'>('style')
  const animationRef = useRef<number | null>(null)
  const lastFrameRef = useRef<number>(0)
  const currentFrameRef = useRef(0)
  const lastAnimTimeRef = useRef<number>(0)

  useEffect(() => { currentFrameRef.current = trajState.currentFrame }, [trajState.currentFrame])

  // Apply style and color to viewer
  const applyStyle = (newStyle: string, newColor: string) => {
    if (!viewerRef.current) return
    
    const viewer = viewerRef.current
    
    // Chain-specific colors for colorful display
    const chainColors: Record<string, string> = {
      A: '#FF6B6B', // Red
      B: '#4ECDC4', // Teal
      C: '#45B7D1', // Blue
      D: '#96CEB4', // Green
      E: '#FFEAA7', // Yellow
      F: '#DDA0DD', // Plum
    }

    if (newColor === 'chainHetatm') {
      // Color by chain with different colors
      viewer.setStyle({}, {
        [newStyle]: {
          colorscheme: {
            prop: 'chain',
            palette: chainColors
          }
        }
      })
    } else if (newColor === 'chain') {
      // Standard chain coloring
      viewer.setStyle({}, {
        [newStyle]: {
          colorscheme: {
            prop: 'chain'
          }
        }
      })
    } else {
      // Standard coloring schemes
      viewer.setStyle({}, {
        [newStyle]: {
          colorscheme: newColor
        }
      })
    }
    
    viewer.render()
  }

  useEffect(() => {
    if (!containerRef.current || !window.$3Dmol) return

    const viewer = window.$3Dmol.createViewer(containerRef.current, {
      backgroundColor: '#1a1a2e', // Dark background for better visibility
    })

    viewerRef.current = viewer

    // Add the protein structure
    viewer.addModel(SAMPLE_PROTEIN, 'pdb')
    
    // Apply initial style with chain coloring
    applyStyle('cartoon', 'chain')
    
    viewer.zoomTo()
    viewer.render()
    setLoaded(true)

    return () => {
      viewer.clear()
    }
  }, [])

  // Load trajectory from URL params
  useEffect(() => {
    const params = new URLSearchParams(window.location.search)
    const trajPath = params.get('trajectory')
    const topoPath = params.get('topology')
    const rmsdPath = params.get('rmsd')

    if (trajPath && containerRef.current && window.$3Dmol) {
      setTrajectoryPath(trajPath)
      setTopologyPath(topoPath)
      setShowTrajectoryControls(true)

      const viewer = viewerRef.current
      if (!viewer) return

      try {
        viewer.addTrajectory(trajPath, topoPath || SAMPLE_PROTEIN)
        const numFrames = viewer.numFrames ? viewer.numFrames() : 100
        
        let rmsdHistory: number[] = []
        if (rmsdPath) {
          try {
            fetch(rmsdPath)
              .then((r) => r.text())
              .then((text) => {
                rmsdHistory = text
                  .split('\n')
                  .filter((l) => l.trim() && !isNaN(parseFloat(l)))
                  .map((l) => parseFloat(l))
                setTrajState((prev) => ({
                  ...prev,
                  rmsdHistory: rmsdHistory.slice(0, numFrames),
                }))
              })
              .catch(() => {})
          } catch {}
        }

        setTrajState((prev) => ({
          ...prev,
          totalFrames: numFrames,
          rmsdHistory: rmsdHistory.slice(0, numFrames),
        }))
        viewer.setFrame(0)
        viewer.render()
        lastFrameRef.current = 0
      } catch (err) {
        console.error('Failed to load trajectory:', err)
      }
    }
  }, [])

  // Update style when selection changes
  useEffect(() => {
    if (viewerRef.current) {
      applyStyle(style, colorScheme)
    }
  }, [style, colorScheme])

  const handleRotate = () => {
    if (!viewerRef.current) return
    const viewer = viewerRef.current
    let rotating = true
    const rotate = () => {
      if (!rotating || !viewerRef.current) return
      viewer.rotate(2, { x: 0, y: 1, z: 0 })
      viewer.render()
      requestAnimationFrame(rotate)
    }
    rotate()
    setTimeout(() => { rotating = false }, 5000)
  }

  const handleZoomFit = () => {
    if (!viewerRef.current) return
    viewerRef.current.zoomTo()
    viewerRef.current.render()
  }

  const handleReset = () => {
    if (!viewerRef.current) return
    viewerRef.current.reset()
    viewerRef.current.zoomTo()
    applyStyle(style, colorScheme)
  }

  const handleScreenshot = () => {
    if (!viewerRef.current) return
    const png = viewerRef.current.png()
    const link = document.createElement('a')
    link.download = `molecule_${Date.now()}.png`
    link.href = png
    link.click()
  }

  const handlePlayPause = () => {
    if (!viewerRef.current) return

    if (trajState.isPlaying) {
      if (animationRef.current) {
        cancelAnimationFrame(animationRef.current)
        animationRef.current = null
      }
      setTrajState((prev) => ({ ...prev, isPlaying: false }))
    } else {
      lastAnimTimeRef.current = 0
      const animate = (timestamp: number) => {
        if (!viewerRef.current) return
        if (lastAnimTimeRef.current === 0) {
          lastAnimTimeRef.current = timestamp
        }
        const elapsed = timestamp - lastAnimTimeRef.current
        const frameInterval = 1000 / trajState.playbackSpeed
        if (elapsed >= frameInterval) {
          lastAnimTimeRef.current = timestamp - (elapsed % frameInterval)
          const nextFrame = (currentFrameRef.current + 1) % trajState.totalFrames
          viewerRef.current.setFrame(nextFrame)
          viewerRef.current.render()
          currentFrameRef.current = nextFrame
          setTrajState((prev) => ({ ...prev, currentFrame: nextFrame }))
        }
        animationRef.current = requestAnimationFrame(animate)
      }
      animationRef.current = requestAnimationFrame(animate)
      setTrajState((prev) => ({ ...prev, isPlaying: true }))
    }
  }

  const handleFrameChange = (frame: number) => {
    if (!viewerRef.current) return
    viewerRef.current.setFrame(frame)
    viewerRef.current.render()
    lastFrameRef.current = frame
    const time_ps = Math.round(frame * 10)
    const rmsd = trajState.rmsdHistory[frame] ?? null
    setTrajState((prev) => ({
      ...prev,
      currentFrame: frame,
      time_ps,
      rmsdHistory: prev.rmsdHistory.map((r, i) => (i === frame && rmsd !== null ? rmsd : r)),
    }))
  }

  const handleSpeedChange = (speed: number) => {
    setTrajState((prev) => ({ ...prev, playbackSpeed: speed }))
  }

  const handleCenterBindingSite = () => {
    if (!viewerRef.current || !bindingSite) return
    const viewer = viewerRef.current
    viewer.centerOn(
      { x: bindingSite.center[0], y: bindingSite.center[1], z: bindingSite.center[2] },
      bindingSite.radius
    )
    viewer.render()
  }

  const handleToggleSurface = () => {
    if (!viewerRef.current) return
    const viewer = viewerRef.current
    if (showSurface) {
      viewer.removeSurface('surface')
    } else {
      viewer.addSurface('surface', { opacity: 0.5, color: 'white' })
    }
    viewer.render()
    setShowSurface(!showSurface)
  }

  const handlePDBUpload = async (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0]
    if (!file || !viewerRef.current) return

    try {
      const content = await file.text()
      const viewer = viewerRef.current
      
      viewer.removeAllModels()
      viewer.addModel(content, 'pdb')
      applyStyle(style, colorScheme)
      viewer.zoomTo()
      viewer.render()
      
      setLoadedFileName(file.name)
      setLoaded(true)
    } catch (err) {
      console.error('Failed to load PDB file:', err)
    }
  }

  return (
    <div className="p-6">
      <div className="mb-6">
        <div className="flex items-center gap-3">
          <h1 className="text-2xl font-bold text-text-primary">3D Molecular Viewer</h1>
          <div className="relative group">
            <span className="text-gray-400 cursor-help">ℹ️</span>
            <div className="absolute bottom-full left-1/2 -translate-x-1/2 mb-2 px-3 py-2 bg-gray-900 text-white text-xs rounded-lg w-72 opacity-0 invisible group-hover:opacity-100 group-hover:visible transition-all z-50 shadow-lg">
              Visualize protein structures in 3D. Drag to rotate, scroll to zoom, right-drag to pan.
              <div className="absolute top-full left-1/2 -translate-x-1/2 border-4 border-transparent border-t-gray-900"></div>
            </div>
          </div>
        </div>
        <p className="text-text-secondary mt-1">
          {showTrajectoryControls
            ? `📽️ Trajectory Playback | Frame ${trajState.currentFrame + 1} of ${trajState.totalFrames || 1}${trajState.rmsdHistory[trajState.currentFrame] !== undefined ? ` | RMSD: ${trajState.rmsdHistory[trajState.currentFrame].toFixed(3)} Å` : ''}`
            : '🧬 Interactive 3D visualization of molecular structures'}
        </p>
      </div>

      <div className="grid grid-cols-1 lg:grid-cols-4 gap-6">
        {/* Viewer */}
        <div className="lg:col-span-3">
          <Card padding="none" className="overflow-hidden">
            {/* Toolbar */}
            <div className="flex flex-wrap items-center gap-2 p-3 bg-surface-secondary border-b border-border-light">
              <span className="text-xs text-text-tertiary mr-2">Controls:</span>
              
              {/* Upload PDB Button */}
              <label className="flex items-center gap-1 px-2 py-1 bg-blue-600 hover:bg-blue-700 text-white rounded text-xs cursor-pointer">
                📤 Upload PDB
                <input type="file" accept=".pdb,.ent" className="hidden" onChange={handlePDBUpload} />
              </label>
              
              {[
                { icon: '🔄', label: 'Auto Rotate', handler: handleRotate },
                { icon: '⊕', label: 'Zoom Fit', handler: handleZoomFit },
                { icon: '↺', label: 'Reset View', handler: handleReset },
                ...(showTrajectoryControls ? [{ icon: '🎬', label: 'Surface', handler: handleToggleSurface }] : []),
                ...(bindingSite ? [{ icon: '🎯', label: 'Focus Site', handler: handleCenterBindingSite }] : []),
                { icon: '📷', label: 'Screenshot', handler: handleScreenshot },
              ].map((btn) => (
                <Button
                  key={btn.label}
                  variant="outline"
                  size="sm"
                  onClick={btn.handler}
                >
                  {btn.icon} {btn.label}
                </Button>
              ))}
            </div>

            {showTrajectoryControls && (
              <div className="flex items-center gap-4 p-3 bg-surface-secondary border-b border-border-light flex-wrap">
                <Button variant="outline" size="sm" onClick={handlePlayPause}>
                  {trajState.isPlaying ? '⏸ Pause' : '▶ Play'}
                </Button>

                <div className="flex items-center gap-2 flex-1 min-w-[200px]">
                  <span className="text-xs text-text-secondary">Frame:</span>
                  <input
                    type="range"
                    min={0}
                    max={trajState.totalFrames > 0 ? trajState.totalFrames - 1 : 0}
                    value={trajState.currentFrame}
                    onChange={(e) => handleFrameChange(parseInt(e.target.value))}
                    className="flex-1"
                  />
                  <span className="text-xs text-text-secondary">
                    {trajState.currentFrame + 1} / {trajState.totalFrames || 1}
                  </span>
                </div>

                {trajState.rmsdHistory.length > 0 && (
                  <div className="flex items-center gap-2 px-3 py-1 bg-primary-50 rounded text-xs">
                    <span className="text-text-secondary">RMSD:</span>
                    <span className="font-mono font-bold text-primary">
                      {trajState.rmsdHistory[trajState.currentFrame]?.toFixed(3) ?? '--'} Å
                    </span>
                  </div>
                )}

                <div className="flex items-center gap-2">
                  <span className="text-xs text-text-secondary">Speed:</span>
                  <select
                    value={trajState.playbackSpeed}
                    onChange={(e) => handleSpeedChange(parseInt(e.target.value))}
                    className="text-xs border rounded px-1"
                  >
                    <option value={5}>5 fps</option>
                    <option value={10}>10 fps</option>
                    <option value={15}>15 fps</option>
                    <option value={30}>30 fps</option>
                  </select>
                </div>
              </div>
            )}

            {/* 3D Viewer */}
            <div
              ref={containerRef}
              className="w-full h-[500px]"
              style={{ position: 'relative', background: '#1a1a2e' }}
            />

            {/* Status */}
            <div className="flex items-center justify-between p-3 bg-surface-secondary border-t border-border-light text-xs text-text-secondary">
              <span>
                {loaded
                  ? showTrajectoryControls
                    ? `Trajectory: ${trajState.totalFrames} frames | Time: ${trajState.time_ps} ps | ${trajState.currentFrame + 1}/${trajState.totalFrames || 1}`
                    : loadedFileName
                      ? `🧬 Loaded: ${loadedFileName}`
                      : '🧬 Sample structure - 2 chains (A & B)'
                  : 'Loading...'}
              </span>
              <span>Use mouse to rotate/zoom | Scroll to zoom</span>
            </div>
          </Card>
        </div>

        {/* Style Panel */}
        <Card>
          <h3 className="font-bold text-text-primary mb-4">Display Style</h3>

          <div className="space-y-2 mb-6">
            {VIEW_STYLES.map((s) => (
              <Button
                key={s.id}
                variant={style === s.id ? 'primary' : 'outline'}
                size="sm"
                className="w-full justify-start"
                onClick={() => setStyle(s.id)}
              >
                {s.label}
              </Button>
            ))}
          </div>

          <h3 className="font-bold text-text-primary mb-4">Color Scheme</h3>
          <div className="space-y-2">
            {COLOR_SCHEMES.map((c) => (
              <Button
                key={c.id}
                variant={colorScheme === c.id ? 'primary' : 'outline'}
                size="sm"
                className="w-full justify-start"
                onClick={() => setColorScheme(c.id)}
              >
                {c.label}
              </Button>
            ))}
          </div>

          {showTrajectoryControls && trajState.rmsdHistory.length > 0 && (
            <div className="mt-6">
              <div className="flex items-center gap-2 mb-3">
                <h3 className="font-bold text-text-primary">RMSD Analysis</h3>
                <div className="relative group">
                  <span className="text-gray-400 cursor-help text-xs">ℹ️</span>
                  <div className="absolute bottom-full left-1/2 -translate-x-1/2 mb-2 px-3 py-2 bg-gray-900 text-white text-xs rounded-lg w-64 opacity-0 invisible group-hover:opacity-100 group-hover:visible transition-all z-50 shadow-lg">
                    RMSD (Root Mean Square Deviation) measures how much the structure moves from the reference frame. Lower values = more stable.
                    <div className="absolute top-full left-1/2 -translate-x-1/2 border-4 border-transparent border-t-gray-900"></div>
                  </div>
                </div>
              </div>
              <div className="p-3 bg-surface-secondary rounded-lg">
                <div className="text-xs space-y-2">
                  <div className="flex justify-between">
                    <span className="text-text-secondary">Current:</span>
                    <span className="font-mono font-bold text-primary">
                      {trajState.rmsdHistory[trajState.currentFrame]?.toFixed(3) ?? '--'} Å
                    </span>
                  </div>
                  <div className="flex justify-between">
                    <span className="text-text-secondary">Min:</span>
                    <span className="font-mono text-green-600">
                      {Math.min(...trajState.rmsdHistory).toFixed(3)} Å
                    </span>
                  </div>
                  <div className="flex justify-between">
                    <span className="text-text-secondary">Max:</span>
                    <span className="font-mono text-red-600">
                      {Math.max(...trajState.rmsdHistory).toFixed(3)} Å
                    </span>
                  </div>
                  <div className="flex justify-between">
                    <span className="text-text-secondary">Avg:</span>
                    <span className="font-mono">
                      {(
                        trajState.rmsdHistory.reduce((a, b) => a + b, 0) /
                        trajState.rmsdHistory.length
                      ).toFixed(3)}{' '}
                      Å
                    </span>
                  </div>
                </div>
                <div className="mt-3 h-16 bg-gray-700 rounded overflow-hidden">
                  <svg viewBox="0 0 100 40" className="w-full h-full" preserveAspectRatio="none">
                    <polyline
                      fill="none"
                      stroke="#4ade80"
                      strokeWidth="1"
                      points={trajState.rmsdHistory
                        .map((v, i) => {
                          const maxRmsd = Math.max(...trajState.rmsdHistory, 0.001)
                          const x = (trajState.rmsdHistory.length > 1 ? i / (trajState.rmsdHistory.length - 1) : 0.5) * 100
                          return `${x},${30 - (v / maxRmsd) * 25}`
                        })
                        .join(' ')}
                    />
                    <line x1="0" y1="30" x2="100" y2="30" stroke="#374151" strokeWidth="0.5" />
                  </svg>
                </div>
                <p className="text-xs text-text-tertiary mt-2">
                  RMSD plot showing structural deviation over time
                </p>
              </div>
            </div>
          )}
        </Card>
      </div>
    </div>
  )
}
