/*

TERMS OF USAGE

This software is Copyright (C) 2008 by Andrew E. Firth, University College
  Cork, Cork, Ireland. All rights reserved.
Updates are Copyright (C) 2014 by Andrew E. Firth, University of Cambridge,
 Cambridge, UK. All rights reserved.

The Software is currently supplied by the Authors free of charge and without
  any guarantee of fitness or usability for any purpose whatsoever. In these
  Terms of Usage, the term "the Software" is considered to include all software
  and data associated with this software and the associated website.

This Software is provided "as is" without warranty of any kind, either express
  or implied, including, without limitation, any warranty of merchantability
  and fitness for a particular purpose. The entire risk as to the results and
  performance of this Software is assumed by you. The author and other parties
  assume no responsibility for the accuracy or application of or errors or
  omissions in this Software. In no event shall the author or any other party
  be liable for any direct, indirect, special, incidental or consequential
  damages arising out of the use or inability to use this Software, even if the
  author or other parties have been advised of the likelihood of such damages
  occurring. The author or any other party shall not be liable for any loss,
  damages or costs, arising out of, but not limited to, lost profits or
  revenue, loss of use of the Software, loss of data or equipment, the costs of
  recovering software, data or equipment, the cost of substitute software or
  data, claims by third parties, or other similar costs.

The author or other parties associated with this Software cannot be held
  responsible or liable for any damage or loss to any equipment, programs or
  data files, time spent or any other expense incured due to the usage of this
  Software.

If you do not agree with these terms do not install this Software on any
  computer system or utilize it. Use of this Software will be taken as an
  agreement to these terms.

Under these conditions, permission is granted to anyone to use this software
  for academic purposes only. The user is prohibited from re-distributing the
  software without permission from the author. The software or any derivatives
  of the software should under no circumstance be re-distributed commercially
  or for profit.

*/

#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <cstring>
using namespace std;
#define maxlength 50000
#define maxpairs 1000

// Function declarations
int readsequence(char seq[128], int gsequences[][2], int seqid);
int calc_nuc(double nuc[][4], double ttt, double nucB1[][4], double nucB2[][4],
	     double nucD[][4]);

int main(int argc, char* argv[])
{
  if (2 != argc) {
    cerr << "Usage: '" << argv[0] << " pairs_file'.\n";
    exit(EXIT_FAILURE);
  }

  //---------------------------------------------------------------------------
  // Variable declarations.

  int g, i, j, k, firstseq, nseq, npairs, alnlength, codonlen, maxiter;
  int np, aa2aa[64], gsequences[maxlength][2], niter, ndiff, test;
  int jj[3][2], jcid[2], k0, k1, k2, kk[3], kcid, flag[maxpairs];
  int codoncount[maxpairs], nmuts[maxpairs], ntransitions[maxpairs];
  double ttt, tttmin, tttmax, marprob, prob, ettt[2], enmuts[2], expmuts;
  double tttx[maxpairs], expnmuts[maxpairs], expntransitions[maxpairs];
  double plotdata[int(maxlength/3)+1][3], tempA[4][4], testA[4][4], tfit;
  double nucA[4][4], nucB1[4][4], nucB2[4][4], nucD[4][4], nuc[4][4];
  double sumsqs, mean;
  char aa2codon[64], seq[128], id[2][128];

  int verbose = 0;

  //---------------------------------------------------------------------------
  // Setup parameters.
  tttmin = 0.0;
  tttmax = 10.0;
  tfit = 0.2;
  maxiter = 20;

  if (verbose) {
    cout << "\n";
  }

  for (i = 0; i < int(maxlength/3)+1; ++i) {
    plotdata[i][0] = plotdata[i][1] = plotdata[i][2] = 0.;
  }

  //---------------------------------------------------------------------------

  // Codon to amino acid table; indexed 0 to 64 as follows:
  /*
  F UUU 00   S UCU 04   Y UAU 08   C UGU 12
  F UUC 01   S UCC 05   Y UAC 09   C UGC 13
  L UUA 02   S UCA 06   Z UAA 10   Z UGA 14
  L UUG 03   S UCG 07   Z UAG 11   W UGG 15
  L CUU 16   P CCU 20   H CAU 24   R CGU 28
  L CUC 17   P CCC 21   H CAC 25   R CGC 29
  L CUA 18   P CCA 22   Q CAA 26   R CGA 30
  L CUG 19   P CCG 23   Q CAG 27   R CGG 31
  I AUU 32   T ACU 36   N AAU 40   S AGU 44
  I AUC 33   T ACC 37   N AAC 41   S AGC 45
  I AUA 34   T ACA 38   K AAA 42   R AGA 46
  M AUG 35   T ACG 39   K AAG 43   R AGG 47
  V GUU 48   A GCU 52   G GGU 60   D GAU 56
  V GUC 49   A GCC 53   G GGC 61   D GAC 57
  V GUA 50   A GCA 54   G GGA 62   E GAA 58
  V GUG 51   A GCG 55   G GGG 63   E GAG 59
  */

  aa2codon[0] = 'F';      aa2codon[32] = 'I';
  aa2codon[1] = 'F';	  aa2codon[33] = 'I';
  aa2codon[2] = 'L';	  aa2codon[34] = 'I';
  aa2codon[3] = 'L';	  aa2codon[35] = 'M';
  aa2codon[4] = 'S';	  aa2codon[36] = 'T';
  aa2codon[5] = 'S';	  aa2codon[37] = 'T';
  aa2codon[6] = 'S';	  aa2codon[38] = 'T';
  aa2codon[7] = 'S';	  aa2codon[39] = 'T';
  aa2codon[8] = 'Y';	  aa2codon[40] = 'N';
  aa2codon[9] = 'Y';	  aa2codon[41] = 'N';
  aa2codon[10] = 'Z';	  aa2codon[42] = 'K';
  aa2codon[11] = 'Z';	  aa2codon[43] = 'K';
  aa2codon[12] = 'C';	  aa2codon[44] = 'S';
  aa2codon[13] = 'C';	  aa2codon[45] = 'S';
  aa2codon[14] = 'Z';	  aa2codon[46] = 'R';
  aa2codon[15] = 'W';	  aa2codon[47] = 'R';
  aa2codon[16] = 'L';	  aa2codon[48] = 'V';
  aa2codon[17] = 'L';	  aa2codon[49] = 'V';
  aa2codon[18] = 'L';	  aa2codon[50] = 'V';
  aa2codon[19] = 'L';	  aa2codon[51] = 'V';
  aa2codon[20] = 'P';	  aa2codon[52] = 'A';
  aa2codon[21] = 'P';	  aa2codon[53] = 'A';
  aa2codon[22] = 'P';	  aa2codon[54] = 'A';
  aa2codon[23] = 'P';	  aa2codon[55] = 'A';
  aa2codon[24] = 'H';	  aa2codon[56] = 'D';
  aa2codon[25] = 'H';	  aa2codon[57] = 'D';
  aa2codon[26] = 'Q';	  aa2codon[58] = 'E';
  aa2codon[27] = 'Q';	  aa2codon[59] = 'E';
  aa2codon[28] = 'R';	  aa2codon[60] = 'G';
  aa2codon[29] = 'R';	  aa2codon[61] = 'G';
  aa2codon[30] = 'R';	  aa2codon[62] = 'G';
  aa2codon[31] = 'R';	  aa2codon[63] = 'G';


  // 4 x 4 diagonalized nucleotide mutation matrix. Order is UCAG. 
  //   row i -> column k.  Four matrices: A, B, D, B^{-1}, where
  //   A = B.D.B^{-1} and D is diagonal.  A itself isn't actually used.

  nucA[0][0] = -1.0; nucA[0][1] =  0.6; nucA[0][2] =  0.2; nucA[0][3] =  0.2;
  nucA[1][0] =  0.6; nucA[1][1] = -1.0; nucA[1][2] =  0.2; nucA[1][3] =  0.2;
  nucA[2][0] =  0.2; nucA[2][1] =  0.2; nucA[2][2] = -1.0; nucA[2][3] =  0.6;
  nucA[3][0] =  0.2; nucA[3][1] =  0.2; nucA[3][2] =  0.6; nucA[3][3] = -1.0;

  nucB1[0][0] = -0.642685; nucB1[0][1] = -0.565194; nucB1[0][2] =  0.0359716; nucB1[0][3] = -0.688247;
  nucB1[1][0] = -0.642685; nucB1[1][1] = -0.565194; nucB1[1][2] = -0.0359716; nucB1[1][3] =  0.688247;
  nucB1[2][0] = -0.642685; nucB1[2][1] =  0.565194; nucB1[2][2] = -0.5;       nucB1[2][3] = -6.48374e-17;
  nucB1[3][0] = -0.642685; nucB1[3][1] =  0.565194; nucB1[3][2] =  0.5;       nucB1[3][3] = -1.10677e-17;

  nucD[0][0] = 0.0; nucD[0][1] =  0.0; nucD[0][2] =  0.0; nucD[0][3] =  0.0;
  nucD[1][0] = 0.0; nucD[1][1] = -0.8; nucD[1][2] =  0.0; nucD[1][3] =  0.0;
  nucD[2][0] = 0.0; nucD[2][1] =  0.0; nucD[2][2] = -1.6; nucD[2][3] =  0.0;
  nucD[3][0] = 0.0; nucD[3][1] =  0.0; nucD[3][2] =  0.0; nucD[3][3] = -1.6;

  nucB2[0][0] = -0.388993;    nucB2[0][1] = -0.388993;    nucB2[0][2] = -0.388993;  nucB2[0][3] = -0.388993;
  nucB2[1][0] = -0.442326;    nucB2[1][1] = -0.442326;    nucB2[1][2] =  0.442326;  nucB2[1][3] =  0.442326;
  nucB2[2][0] =  8.06558e-17; nucB2[2][1] = -8.06558e-17; nucB2[2][2] = -1.0;       nucB2[2][3] =  1.0;
  nucB2[3][0] = -0.726483;    nucB2[3][1] =  0.726483;    nucB2[3][2] = -0.0522655; nucB2[3][3] =  0.0522655;

  // Check B.D.B{^-1}.
  for (i = 0; i < 4; ++i) {
    for (j = 0; j < 4; ++j) {
      tempA[i][j] = nucB1[i][j]*nucD[j][j];
    }
  }
  for (i = 0; i < 4; ++i) {
    for (j = 0; j < 4; ++j) {
      testA[i][j] = 0.;
      for (k = 0; k < 4; ++k) {
	testA[i][j] += tempA[i][k]*nucB2[k][j];
      }
    }
  }
  if (verbose) {
    cout << "Product B D B^{-1}\n";
    for (i = 0; i < 4; ++i) {
      cout << "   ";
      for (j = 0; j < 4; ++j) {
	cout << testA[i][j] << " ";
      }
      cout << "\n";
    }
    cout << "\n";
  }
  for (i = 0; i < 4; ++i) {
    for (k = 0; k < 4; ++k) {
      if (i == k) {
	if (testA[i][k] > 0.) {
	  cerr << "Aborting: Nucleotide matrix value at " << i+1 << ", " 
	       << k+1 << " is positive.\n";
	  exit(EXIT_FAILURE);
	}
      } else {
	if (testA[i][k] < 0.) {
	  cerr << "Aborting: Nucleotide matrix value at " << i+1 << ", " 
	       << k+1 << " is negative.\n";
	  exit(EXIT_FAILURE);
	}
      }
    }
  }

  //---------------------------------------------------------------------------
  // Miscellaneous stuff.

  // Fitting parameters for ttt (evolutionary time).
  if (tttmin < 0.) {
    tttmin = 0.;
  }
  if (tttmax < tttmin) {
    tttmax = tttmin;
  }
  if (tfit <= 0 || tfit > 1.) {
    tfit = 1.;
  }
  if (maxiter < 1) {
    maxiter = 1;
  }

  // Number of sequence pairs in sequences file.
  ifstream sequencesfile(argv[1]);
  if (!sequencesfile) {
    cerr << "Aborting: Can't find pairs file '" << argv[1] << "'.\n";
    exit(EXIT_FAILURE);
  }
  npairs = -1;
  while (sequencesfile.ignore(1000, '\n')) {
    ++npairs;
  }
  sequencesfile.clear();
  sequencesfile.seekg(0);

  //---------------------------------------------------------------------------
  // Calculate conservation scores.

  // Open output data file (one line per sequence pair).
  ofstream outputfile("synplot_data.txt");
  if (!outputfile) {
    cerr << "Failed to open output file 'synplot_data.txt'.\n";
    exit(EXIT_FAILURE);
  }

  firstseq = 1;   // Find alnlength from first sequence read in.

  // Loop through sequence pairs.
  for (np = 0; np < npairs; ++np) {

    // Read sequence pair.
    for (k = 0; k < 2; ++k) {
      if (0 == k) {
	sequencesfile.getline(id[k], 128, ' ');
	cerr << "Processing sequence pair " << id[k] << " + ";
      } else {
	sequencesfile.getline(id[k], 128, '\n');
	cerr << id[k] << "\n";
      }

      // Read sequence pair.
      strcpy(seq, id[k]);
      nseq = readsequence(seq, gsequences, k);
      if (verbose) {
	cout << "Read '" << seq << "', length = " << nseq << ".\n";
      }
      if (firstseq) {
	alnlength = nseq;
	if (0 != alnlength % 3) {
	  cerr << "Aborting: Alignment length not a multiple of 3 nt. "
	       << "You must use a back-translated amino acid alignment.\n";
	  exit(EXIT_FAILURE);
	}
	codonlen = alnlength / 3;
	firstseq = 0;
      }
      if (nseq != alnlength) {
	cerr << "Aborting: Sequences of differing length."
	     << " Input sequences should have alignment gaps (i.e. '-'s)"
	     << " so that they are all the same length.\n";
	exit(EXIT_FAILURE);
      }

      // Check for gaps not in groups of 3 in coding sequences
      for (i = 0; i < codonlen; ++i) {
	if (99 == gsequences[3*i][k] || 99 == gsequences[3*i+1][k] 
	    || 99 == gsequences[3*i+2][k]) {
	  if (99 != gsequences[3*i][k] || 99 != gsequences[3*i+1][k] 
	      || 99 != gsequences[3*i+2][k]) {
	    cerr << "Aborting: Found partially gapped codon - sequence pair "
		 << np+1 << ", sequence " << k+1 << ", alignment coords " 
		 << 3*i+1 << ".." << 3*i+3 << ".\n";
	    cerr << "Gaps must occur in groups of three and in the zero-frame."
		 << "\n";
	    exit(EXIT_FAILURE);
	  }
	}
      }

    } // Close loop for reading in the two sequences in a sequence pair.

    // Calculate observed stats for the sequence pair
    codoncount[np] = 0;   // number of synonymous codon sites (no ambig nt)
                          //   (x 3 to get number of nt in synonymous sites)
    nmuts[np] = 0;        // number of nt mutations within synonymous sites
    ntransitions[np] = 0; // number of transitions within synonymous sites
                          //   (ntransversions = nmuts - ntransitions)
    for (i = 0; i < codonlen; ++i) {
      for (k = 0; k < 2; ++k) {
	for (j = 0; j < 3; ++j) {
	  jj[j][k] = gsequences[3*i+j][k];
	}
	jcid[k] = 16 * jj[0][k] + 4 * jj[1][k] + jj[2][k];
      }
      if (jcid[0] < 64 && jcid[1] < 64) {  // No gap or ambig nt
	if (aa2codon[jcid[0]] == aa2codon[jcid[1]]) { // Synonymous site
	  codoncount[np] += 1;
	  for (j = 0; j < 3; ++j) {
	    if (jj[j][0] != jj[j][1]) {
	      nmuts[np] += 1;
	      plotdata[i][0] += 1; // Observed syn muts summed over aln
	      if ((0 == jj[j][0] && 1 == jj[j][1])
		  ||(1 == jj[j][0] && 0 == jj[j][1])
		  ||(2 == jj[j][0] && 3 == jj[j][1])
		  ||(3 == jj[j][0] && 2 == jj[j][1])) {
		ntransitions[np] += 1;
	      }
	    }
	  }
	}
      }
    }
    if (verbose) {
      cerr << "Sequence pair " << np+1 << ": there are " << codoncount[np] 
	   << " synonymous codon sites, and " << nmuts[np] 
	   << " point mutations in synonymous sites, of which " 
	   << ntransitions[np] << " are transitions.\n";
    }

    // Find ttt value giving closest match between expected # muts and observed
    //   # muts (at synonymous sites).  Use bisection method -> since nmuts 
    //   increases monotonically with ttt.  Stop when abs(exp # muts - obs #
    //   muts) < tfit.  Exit if number of iterations exceeds maxiter or if obs
    //   # muts is outside range allowed by tttmin--tttmax.
    flag[np] = 0;         // flag = 1 means problem with ttt fit convergence

    ettt[0] = tttmin; 
    ettt[1] = tttmax;

    // Calculate exp#muts at syn sites for ettt[].
    for (g = 0; g < 2; ++g) {
      test = calc_nuc(nuc, ettt[g], nucB1, nucB2, nucD);
      enmuts[g] = 0.;
      for (i = 0; i < codonlen; ++i) {
	for (k = 0; k < 2; ++k) {
	  for (j = 0; j < 3; ++j) {
	    jj[j][k] = gsequences[3*i+j][k];
	  }
	  jcid[k] = 16 * jj[0][k] + 4 * jj[1][k] + jj[2][k];
	}
	if (jcid[0] < 64 && jcid[1] < 64) {  // No gap or ambig nt
	  if (aa2codon[jcid[0]] == aa2codon[jcid[1]]) { // Synonymous site
	    for (k = 0; k < 2; ++k) { // Fwd and rev directions
	      // First calculate the marginal prob, i.e. the prob that jcid[k]
	      //   mutates to any synonymous codon given nuc[][].  Since we
	      //   know a priori that the site is synonymous, we need to 
	      //   scale probs for individual jcid[k] -> kcid (syn) by this
	      //   marginal prob.
	      marprob = 0.;
	      for (k0 = 0; k0 < 4; ++k0) {
		for (k1 = 0; k1 < 4; ++k1) {
		  for (k2 = 0; k2 < 4; ++k2) {
		    kcid = 16 * k0 + 4 * k1 + k2; // All possible codons
		    if (aa2codon[kcid] == aa2codon[jcid[0]]) { // Synonymous
		      // Probability of jcid[k] mutating to kcid given nuc[][]
		      marprob += nuc[jj[0][k]][k0] * nuc[jj[1][k]][k1] 
			* nuc[jj[2][k]][k2];
		    }
		  }
		}
	      }
	      // Now calculate the expected number of point muts
	      for (k0 = 0; k0 < 4; ++k0) {
		for (k1 = 0; k1 < 4; ++k1) {
		  for (k2 = 0; k2 < 4; ++k2) {
		    kcid = 16 * k0 + 4 * k1 + k2; // All possible codons
		    if (aa2codon[kcid] == aa2codon[jcid[0]]) { // Synonymous
		      // Number of point differences between kcid and jcid[k] 
		      ndiff = 0;
		      if (k0 != jj[0][k]) {ndiff += 1;}
		      if (k1 != jj[1][k]) {ndiff += 1;}
		      if (k2 != jj[2][k]) {ndiff += 1;}
		      // Probability of jcid[k] mutating to kcid given nuc[][]
		      prob = nuc[jj[0][k]][k0] * nuc[jj[1][k]][k1] 
			* nuc[jj[2][k]][k2] / marprob;
		      // Expected number of nt muts for this codon substitution
		      //   (* 0.5 since counting fwd and rev directions)
		      enmuts[g] += 0.5 * float(ndiff) * prob;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }

    // Check tttmin/tttmax range is sufficient.
    if ((enmuts[0] - tfit) > float(nmuts[np])) {
      cerr << "Fitting number of mutations (" << nmuts[np] << ") outside that "
	   << "allowed by given t range.  Sequence pair " << np + 1 << ".\n";
      tttx[np] = tttmin;
      flag[np] = 1;
    } else if (enmuts[1] < float(nmuts[np])) {
      cerr << "Fitting number of mutations (" << nmuts[np] << ") outside that "
	   << "allowed by given t range.  Sequence pair " << np + 1 << ".\n";
      tttx[np] = tttmax;
      flag[np] = 1;
    } else {
	
      // Iterate on ttt.
      niter = 0;
      while (abs(enmuts[0] - float(nmuts[np])) > tfit) {
	
	niter++;
	if (niter > maxiter) {
	  cerr << "Iteration on t failed to converge in " << maxiter 
	       << " iterations.  Sequence pair " << np + 1 << ".\n";
	  flag[np] = 1;
	  break;
	}
	
	// Bisector.
	ttt = 0.5 * (ettt[0] + ettt[1]);
	  
	// Calculate expected number of mutations at bisector.
	test = calc_nuc(nuc, ttt, nucB1, nucB2, nucD);
	expmuts = 0.;
	for (i = 0; i < codonlen; ++i) {
	  for (k = 0; k < 2; ++k) {
	    for (j = 0; j < 3; ++j) {
	      jj[j][k] = gsequences[3*i+j][k];
	    }
	    jcid[k] = 16 * jj[0][k] + 4 * jj[1][k] + jj[2][k];
	  }
	  if (jcid[0] < 64 && jcid[1] < 64) {  // No gap or ambig nt
	    if (aa2codon[jcid[0]] == aa2codon[jcid[1]]) { // Synonymous site
	      for (k = 0; k < 2; ++k) { // Fwd and rev directions
		// First calculate the marginal prob.
		marprob = 0.;
		for (k0 = 0; k0 < 4; ++k0) {
		  for (k1 = 0; k1 < 4; ++k1) {
		    for (k2 = 0; k2 < 4; ++k2) {
		      kcid = 16 * k0 + 4 * k1 + k2; // All possible codons
		      if (aa2codon[kcid] == aa2codon[jcid[0]]) { // Synonymous
			// Prob of jcid[k] mutating to kcid given nuc[][]
			marprob += nuc[jj[0][k]][k0] * nuc[jj[1][k]][k1] 
			  * nuc[jj[2][k]][k2];
		      }
		    }
		  }
		}
		// Now calculate the expected number of point muts
		for (k0 = 0; k0 < 4; ++k0) {
		  for (k1 = 0; k1 < 4; ++k1) {
		    for (k2 = 0; k2 < 4; ++k2) {
		      kcid = 16 * k0 + 4 * k1 + k2; // All possible codons
		      if (aa2codon[kcid] == aa2codon[jcid[0]]) { // Synonymous
			// Number of point differences between kcid and jcid[k] 
			ndiff = 0;
			if (k0 != jj[0][k]) {ndiff += 1;}
			if (k1 != jj[1][k]) {ndiff += 1;}
			if (k2 != jj[2][k]) {ndiff += 1;}
			// Prob of jcid[k] mutating to kcid given nuc[][]
			prob = nuc[jj[0][k]][k0] * nuc[jj[1][k]][k1] 
			  * nuc[jj[2][k]][k2] / marprob;
			// Expected number of nt muts for this substitution
			//   (* 0.5 since counting fwd and rev directions)
			expmuts += 0.5 * float(ndiff) * prob;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
	  
	// Reassign ettt[] and enmuts[].
	if (float(nmuts[np]) > expmuts) {
	  ettt[0] = ttt;
	  enmuts[0] = expmuts;
	} else {
	  ettt[1] = ttt;
	  enmuts[1] = expmuts;
	}
	
      }

      tttx[np] = ettt[0];

    } // close if/else for nmuts[np] outside/within allowed ttt range

    // Calculate nucleotide matrix for best ttt.
    test = calc_nuc(nuc, tttx[np], nucB1, nucB2, nucD);

    // Calculate expected nmuts and ntransitions for best ttt value.
    expnmuts[np] = 0.;
    expntransitions[np] = 0.;
    for (i = 0; i < codonlen; ++i) {
      // These two params summed for a given codon site in a given sequence 
      //   pair.
      sumsqs = 0.;  // Sum of squares, for calculating stddev.
      mean = 0.;    // Mean, for calculating stddev.
      for (k = 0; k < 2; ++k) {
	for (j = 0; j < 3; ++j) {
	  jj[j][k] = gsequences[3*i+j][k];
	}
	jcid[k] = 16 * jj[0][k] + 4 * jj[1][k] + jj[2][k];
      }
      if (jcid[0] < 64 && jcid[1] < 64) {  // No gap or ambig nt
	if (aa2codon[jcid[0]] == aa2codon[jcid[1]]) { // Synonymous site
	  for (k = 0; k < 2; ++k) { // Fwd and rev directions
	    // First calculate the marginal prob.
	    marprob = 0.;
	    for (k0 = 0; k0 < 4; ++k0) {
	      for (k1 = 0; k1 < 4; ++k1) {
		for (k2 = 0; k2 < 4; ++k2) {
		  kcid = 16 * k0 + 4 * k1 + k2; // All possible codons
		  if (aa2codon[kcid] == aa2codon[jcid[0]]) { // Synonymous
		    // Prob of jcid[k] mutating to kcid given nuc[][]
		    marprob += nuc[jj[0][k]][k0] * nuc[jj[1][k]][k1] 
		      * nuc[jj[2][k]][k2];
		  }
		}
	      }
	    }
	    // Now calculate the expected number of point muts
	    for (k0 = 0; k0 < 4; ++k0) {
	      for (k1 = 0; k1 < 4; ++k1) {
		for (k2 = 0; k2 < 4; ++k2) {
		  kcid = 16 * k0 + 4 * k1 + k2; // All possible codons
		  if (aa2codon[kcid] == aa2codon[jcid[0]]) { // Synonymous
		    // Prob of jcid[k] mutating to kcid given nuc[][]
		    prob = nuc[jj[0][k]][k0] * nuc[jj[1][k]][k1] 
		      * nuc[jj[2][k]][k2] / marprob;
		    // Sum up point differences between kcid and jcid[k] 
		    kk[0] = k0; kk[1] = k1; kk[2] = k2;
		    ndiff = 0;
		    for (j = 0; j < 3; ++j) {
		      if (kk[j] != jj[j][k]) {
			ndiff += 1;
			if ((0 == kk[j] && 1 == jj[j][k])
			    ||(1 == kk[j] && 0 == jj[j][k])
			    ||(2 == kk[j] && 3 == jj[j][k])
			    ||(3 == kk[j] && 2 == jj[j][k])) {
			  expntransitions[np] += 0.5 * prob;
			}
		      }
		    }
		    // * 0.5 since counting fwd and rev directions
		    // Expected syn muts summed over aln pair
		    expnmuts[np] += 0.5 * prob * float(ndiff);
		    // Expected syn muts summed over aln col
		    plotdata[i][1] += 0.5 * prob * float(ndiff);
		    // Mean and sum of squares (given codon, given pair, i.e.
		    //   over which Sum(0.5 * prob) = 1).
		    mean += 0.5 * prob * float(ndiff);
		    sumsqs += 0.5 * prob * float(ndiff) * float(ndiff);
		  }
		}
	      }
	    }
	  }
	}
      }
      // Summed over a large column and sliding window the distribution
      //   should eventually be more-or-less approximately normal (by CLT?),
      //   allowing you to use it to calculate p-values in the plot.
      // stddev = sqrt(E(X^2) - E(X)^2) [stddev for the aln col]
      // variance = (E(X^2) - E(X)^2)   [variance for the aln col]
      plotdata[i][2] += (sumsqs - (mean * mean));  // VARIANCE
    }
    if (verbose) {
      cerr << "Sequence pair " << np+1 << ": there are " << expnmuts[np] 
	   << " expected point mutations in synonymous sites, of which " 
	   << expntransitions[np] << " are transitions.\n";
    }

    // Write out stats for each sequence pair.
    outputfile << setprecision(3);
    outputfile << setw(20) << id[0] << " ";
    outputfile << setw(20) << id[1] << " ";
    outputfile << setw(6) << codonlen << " ";
    outputfile << setw(6) << codoncount[np] << " ";
    outputfile << setw(6) << nmuts[np] << " ";
    outputfile << setw(6) << ntransitions[np] << " ";
    outputfile << setw(8) << tttx[np] << " ";
    outputfile << setw(8) << expnmuts[np] << " ";
    outputfile << setw(8) << expntransitions[np] << " ";
    outputfile << setw(6) << flag[np] << "\n";
    /*
    Columns (row for each sequence pair).
    1) Total number of codons (alignment coords).
    2) Number of synonymous sites.
    3) Number of point nt substitutions within these synonymous sites.
    4) Number of transitions.
    5) Best fit t value (more-or-less evolutionary time).
    6) Expected number of point nt substitutions given t.
    7) Expected number of transitions given t.
    8) Flag (0 is OK; 1 means the t convergence failed)
    */
        
    if (verbose) {
      cout << "\n";
    }

  } // Close loop through sequence pairs.

  outputfile.close();
  sequencesfile.close(); 

  // Write out codon-by-codon observed and expected number of base
  //   substitutions at synonymous sites.  Numbers summed over all sequence
  //   pairs.
  ofstream plotfile("synplot_plot.txt");
  if (!plotfile) {
    cerr << "Failed to open log file 'synplot_plot.txt'.\n";
    exit(EXIT_FAILURE);
  }
  for (i = 0; i < codonlen; ++i) {
    plotfile << setprecision(6);
    plotfile << setw(8) << i + 1 << " "
	     << setw(13) << 0.5*plotdata[i][0] << " "
	     << setw(13) << 0.5*plotdata[i][1] << " "
	     << setw(13) << sqrt(0.5*plotdata[i][2]) << "\n";
  }
  plotfile.close();
  /*
    Columns (row for each alignment codon).
    1) Codon index (alignment coords).
    2) obs#muts summed over all pairs
    3) exp#muts summed over all pairs
    4) stddev estimated from exp#muts

    Note that the `0.5*` is to correct for crossing each branch of the tree
    with two pairwise comparisons. So this assumes that the input sequence 
    pairs file traces round the outside of one 2-d realization of the
    phylogenetic tree as recommended.
    
  */

  return(0);
}

//-----------------------------------------------------------------------------

int readsequence(char seq[128], int gsequences[][2], int seqid)

{
  int i, k, nlines, nseq;
  char seqline[maxlength], inputseq[maxlength];

  ifstream sequence(seq);
  if (!sequence) {
    cerr << "Aborting: Can't find sequence file '" << seq << "'.\n";
    exit(EXIT_FAILURE);
  }
  nlines = -1;
  while (sequence.ignore(maxlength, '\n')) {
    ++nlines;
  }
  sequence.clear();
  sequence.seekg(0);
  sequence.ignore(1000, '\n');
  for (i = 0; i < nlines; ++i) {
    sequence.getline(seqline, maxlength, '\n');
    if (0 == i) {
      strcpy(inputseq, seqline);
    } else {
      strcat(inputseq, seqline);
    }
  }
  sequence.close();
  for (i = 0; ; ++i) {
    if (inputseq[i] == 'U' || inputseq[i] == 'u') {
      gsequences[i][seqid] = 0;
    }
    else if (inputseq[i] == 'T' || inputseq[i] == 't') {
      gsequences[i][seqid] = 0;
    }
    else if (inputseq[i] == 'C' || inputseq[i] == 'c') {
      gsequences[i][seqid] = 1;
    }
    else if (inputseq[i] == 'A' || inputseq[i] == 'a') {
      gsequences[i][seqid] = 2;
    }
    else if (inputseq[i] == 'G' || inputseq[i] == 'g') {
      gsequences[i][seqid] = 3;
    }
    else if (inputseq[i] == '-' || inputseq[i] == '.') {
      gsequences[i][seqid] = 99;
    }
    else if (inputseq[i] == '\0') {      
      nseq = i;
      break;
    } else {
      if ((inputseq[i] == 'M' || inputseq[i] == 'm' ||
	   inputseq[i] == 'R' || inputseq[i] == 'r' ||
	   inputseq[i] == 'W' || inputseq[i] == 'w' ||
	   inputseq[i] == 'V' || inputseq[i] == 'v' ||
	   inputseq[i] == 'H' || inputseq[i] == 'h' ||
	   inputseq[i] == 'D' || inputseq[i] == 'd' ||
	   inputseq[i] == 'B' || inputseq[i] == 'b' ||
	   inputseq[i] == 'S' || inputseq[i] == 's' ||
	   inputseq[i] == 'Y' || inputseq[i] == 'y' ||
	   inputseq[i] == 'K' || inputseq[i] == 'k' ||
	   inputseq[i] == 'N' || inputseq[i] == 'n' ||
	   inputseq[i] == 'X' || inputseq[i] == 'x')) {
	gsequences[i][seqid] = 100; // ambiguous nt
	//	cerr << "Warning: Sequence " << seq << ", ambiguous nt '" 
	//		 << inputseq[i] << "', at " << i+1 << ".\n";
      } else {
	cerr << "Aborting: Sequence " << seq << ", unknown nucleotide '" 
		 << inputseq[i] << "', at " << i+1 << ".\n";
	exit(EXIT_FAILURE);
      }
    }
  }
  
  return(nseq);

}


//-----------------------------------------------------------------------------
// Makes 4x4 nucleotide P matrix for a given t.
 
int calc_nuc(double nuc[][4], double ttt, double nucB1[][4], double nucB2[][4],
	     double nucD[][4])

{
  int ni, nj, nk;
  double tempA[4][4];

  for (ni = 0; ni < 4; ++ni) {
    for (nj = 0; nj < 4; ++nj) {
      tempA[ni][nj] = nucB1[ni][nj]*exp(ttt*nucD[nj][nj]);
    }
  }
  for (ni = 0; ni < 4; ++ni) {
    for (nj = 0; nj < 4; ++nj) {
      nuc[ni][nj] = 0.;
      for (nk = 0; nk < 4; ++nk) {
	nuc[ni][nj] += tempA[ni][nk]*nucB2[nk][nj];
      }
      if (nuc[ni][nj] < 0.) { 
	nuc[ni][nj] = 0.; // in case of rounding errors
      }
    }
  }
  return(0);

}

//-----------------------------------------------------------------------------
