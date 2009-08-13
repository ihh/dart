#include "amap/amap_adapter.h"
#include "amap/MultiSequenceDag.h"

#define MIN_PROB 0.01
#define MAX_CONSISTENCY_REPS 5


AMAP_adapter::AMAP_adapter (const FASTA_sequence_database& seq_db, const Dotplot_map& dotplots)
  : seq_db (seq_db),
    num_seqs (seq_db.size()),
    dotplots (dotplots),
    enableVerbose (CTAGGING(3,AMAP_ADAPTER)),
    sequences (seq_db),
    sparseMatrices (num_seqs, SafeVector<SparseMatrix *>(num_seqs, NULL))
{

  // convert dotplots to sparseMatrices
  CTAG(2,AMAP_ADAPTER) << "Converting dotplots to SparseMatrices\n";
  for (int i = 0; i < num_seqs; ++i)
    for (int j = i + 1; j < num_seqs; ++j)
      {
	// Dotplots are in a slightly weird format: x is the horizontal axis => columns
	// and y is the vertical axis => rows
	// It follows that dotplot (x,y) is the yth row and xth column.
	CTAG(2,AMAP_ADAPTER) << "converting for (" << j << "," << i << ")\n";
	if (dotplots.find(i) == dotplots.end()
	    || ((Dotplot_map&)dotplots)[i].find(j) == ((Dotplot_map&)dotplots)[i].end())
	  THROWEXPR ("Can't find post prob matrix (" << i << "," << j << ")");

	VF* posteriorVector = posteriorVector_from_dotplot (((Dotplot_map&)dotplots)[i][j]);
	sparseMatrices[j][i] = new SparseMatrix (sequences.GetSequence(j)->GetLength(), sequences.GetSequence(i)->GetLength(), *posteriorVector);
	sparseMatrices[i][j] = sparseMatrices[j][i]->ComputeTranspose();
      }

  // debugging
  if (CTAGGING(2,AMAP_ADAPTER))
    {
      CL << "Sparse posterior prob. matrices are:\n";
      for (int i = 0; i < num_seqs; ++i)
	for (int j = i + 1; j < num_seqs; ++j)
	  {
	    sparseMatrices[j][i]->Print (CL);
	    sparseMatrices[i][j]->Print (CL);
	  }
    }

}

Stockade AMAP_adapter::get_alignment (float gapFactor /* = 1 */, float edgeWeightThreshold /* = 0 */, int numConsistencyReps /* = 0 */, bool show_intermediates /* = false */, bool output_for_gui /* = false */)
{

  // to do: these should be controlled through parameters to the function
  bool enableEdgeReordering = true;
  bool useTgf = true;

  // control GUI output
  float guiStartWeight = numeric_limits<float>::max();
  int guiStepSize = 1;

  // log posteriors in AMAP format if requested
  if (enableVerbose)
    {
      int numSeqs = sequences.GetNumSequences();
      for (int i = 0; i < numSeqs; i++){
	string seq1name = sequences.GetSequence(i)->GetHeader();
	for (int j = i + 1; j < numSeqs; j++){
	  cout << "Sparse Matrix: " << i << "," << j << endl;
	  cout << "Sequence names: " << seq1name << ", " << sequences.GetSequence(j)->GetHeader() << endl;
	  sparseMatrices[i][j]->Print(cout);
	}
      }
    }

  // perform the consistency transformation if requested
  if (numConsistencyReps > MAX_CONSISTENCY_REPS)
    {
      CTAG(4,AMAP_ADAPTER) << "Reducing requested " << numConsistencyReps << " consistency reps to the max allowed (" << MAX_CONSISTENCY_REPS << ")\n";
      numConsistencyReps = MAX_CONSISTENCY_REPS;
    }
  for (int r = 0; r < numConsistencyReps; r++)
    {
      CTAG(4,AMAP_ADAPTER) << "Performing consistency transformation (" << r << " of " << numConsistencyReps << ")\n";
      // perform consistency transformation to get new posterior matrices
      SafeVector<SafeVector<SparseMatrix *> > newSparseMatrices = DoRelaxation();
      
      // now replace the old posterior matrices
      for (int i = 0; i < num_seqs; i++)
	for (int j = 0; j < num_seqs; j++)
	  {
	    // don't delete the original sparse matrices
	    if (r > 0) delete sparseMatrices[i][j];
	    sparseMatrices[i][j] = newSparseMatrices[i][j];
	  }
    }

  // use AMAP to create an alignment
  CTAG(5,AMAP_ADAPTER) << "Creating a multiple alignment with AMAP.\n";
  bool isNuc = seq_db.alphabet().size() == RNA_alphabet.size() ? true : false;
  if (CTAGGING(5,AMAP_ADAPTER) && output_for_gui) {
    if (isNuc)
      CL << "Creating GUI-format output (nucleotide alphabet).\n";
    else
      CL << "Creating GUI-format output (amino acid alphabet).\n";
  }
  MultiSequenceDag mds (&sequences, false, output_for_gui, guiStartWeight, guiStepSize, isNuc);
  Stockade stockade = mds.AlignDag (sparseMatrices, gapFactor, enableVerbose, enableEdgeReordering, useTgf, edgeWeightThreshold, show_intermediates);

  return stockade;
}
  
VF* AMAP_adapter::posteriorVector_from_dotplot (const Dotplot& dotplot) {

  const int xlen = (int) dotplot.xseq.size();
  const int ylen = (int) dotplot.yseq.size();
  VF *posteriorPtr = new VF ((xlen+1) * (ylen+1), (float) -2e20);
  VF &data = *posteriorPtr;

  // dotplots are in a slightly weird format: x is the horizontal axis => columns
  // and y is the vertical axis => rows
  // It follows that dotplot (x,y) is the yth row and xth column.

  int idx = 0;
  // create leading row of 0's (don't forget that the rows of the SparseMatrix are (xlen+1)-long!!)
  for (int x = 0; x < (int) xlen + 1; ++x)
    data[idx++] = 0.;
  for (int y = 0; y < (int) ylen; ++y)
    {
      // create leading 0
      data[idx++] = 0.;
      for (int x = 0; x < (int) xlen; ++x)
	{
	  // round down to 0
	  if (dotplot (x,y) < MIN_PROB)
	    data[idx++] = 0.;
	  // no probs > 1
	  else if (dotplot (x,y) > 1)
	    data[idx++] = 1.;
	  else
	    data[idx++] = dotplot (x,y);
	}
    }

  // sanity check on dimensions
  if (idx != (xlen+1)*(ylen+1))
    THROWEXPR ("Invalid post prob matrices.");

  return posteriorPtr;
}

/////////////////////////////////////////////////////////////////
// DoRelaxation()
//
// Performs one round of the consistency transformation.  The
// formula used is:
//                     1
//    P'(x[i]-y[j]) = ---  sum   sum P(x[i]-z[k]) P(z[k]-y[j])
//                    |S| z in S  k
//
// where S = {x, y, all other sequences...}
//
/////////////////////////////////////////////////////////////////

SafeVector<SafeVector<SparseMatrix *> > AMAP_adapter::DoRelaxation ()
{
  
  const int numSeqs = sequences.GetNumSequences();

  SafeVector<SafeVector<SparseMatrix *> > newSparseMatrices (numSeqs, SafeVector<SparseMatrix *>(numSeqs, NULL));

  // for every pair of sequences
  for (int i = 0; i < numSeqs; i++){
    for (int j = i+1; j < numSeqs; j++){
      Sequence *seq1 = sequences.GetSequence (i);
      Sequence *seq2 = sequences.GetSequence (j);

      if (enableVerbose)
        cerr << "Relaxing (" << i+1 << ") " << seq1->GetHeader() << " vs. "
             << "(" << j+1 << ") " << seq2->GetHeader() << ": ";

      // get the original posterior matrix
      // (in the form of a matrix coerced into a vector, but with an entry for every possibly residue pairing;
      //  cf the SparseMatrix form)
      VF *posteriorPtr = sparseMatrices[i][j]->GetPosterior(); assert (posteriorPtr);
      VF &posterior = *posteriorPtr;

      const int seq1Length = seq1->GetLength();
      const int seq2Length = seq2->GetLength();

      // oldSums[k] = 1 - gap posterior = "paired posterior"
      //  (the probability that residue k was paired to any residue in the other sequence, ie not paired to a gap)
      VF *oldSumsPtr = new VF(seq1Length + seq2Length + 2,0);
      VF &oldSums = *oldSumsPtr;
      VF *newSumsPtr = new VF(seq1Length + seq2Length + 2,0);
      VF &newSums = *newSumsPtr;

      for (int k = 0, kl = 0; k <= seq1Length; k++) {
	for (int l = 0; l <= seq2Length; l++) {
	  oldSums[k] += posterior[kl];
	  oldSums[seq1Length + 1 + l] += posterior[kl++];
	}
      }

      // contribution from the summation where z = x and z = y
      for (int k = 0; k < (seq1Length+1) * (seq2Length+1); k++) posterior[k] += posterior[k];

      if (enableVerbose)
        cerr << sparseMatrices[i][j]->GetNumCells() << " --> ";

      // contribution from all other sequences
      for (int k = 0; k < numSeqs; k++) if (k != i && k != j){
        Relax (sparseMatrices[i][k], sparseMatrices[k][j], posterior);
      }

      // now renormalization
      for (int k = 0; k < (seq1Length+1) * (seq2Length+1); k++) posterior[k] /= numSeqs;

      for (int k = 0, kl = 0; k <= seq1Length; k++) {
	for (int l = 0; l <= seq2Length; l++) {
	  newSums[k] += posterior[kl];
	  newSums[seq1Length + 1 + l] += posterior[kl++];
	}
      }

      // deal with the gap posteriors
      int gapPostBase = (seq1Length+1) * (seq2Length+1); // offset to beginning of gap posteriors
      for (int k = 0; k < seq1Length + seq2Length + 2; k++) {
	// if paired posterior < POSTERIOR_CUTOFF
	if (oldSums[k] < POSTERIOR_CUTOFF) {
	  if (newSums[k] > 1)
	    cerr << "negative new gap posterior!\n";
	  else {
	    if (enableVerbose)
	      cerr << setprecision(5) << posterior[gapPostBase + k] << "->" << setprecision(5) << 1 - newSums[k] << ", ";
	    posterior[gapPostBase + k] = 1 - newSums[k];
	  }
	}
	else {
	  posterior[gapPostBase + k] *= newSums[k] / oldSums[k];
	  if (enableVerbose && newSums[k] > oldSums[k])
	    cerr << setprecision(5) << newSums[k] / oldSums[k] << ", ";
	}
      }
      
      if (enableVerbose)
	cerr << endl;

      // save the new posterior matrix
      newSparseMatrices[i][j] = new SparseMatrix (seq1->GetLength(), seq2->GetLength(), posterior);
      newSparseMatrices[j][i] = newSparseMatrices[i][j]->ComputeTranspose();

      if (enableVerbose)
        cerr << newSparseMatrices[i][j]->GetNumCells() << " -- ";

      delete posteriorPtr;
      delete oldSumsPtr;
      delete newSumsPtr;

      if (enableVerbose)
        cerr << "done." << endl;
    }
  }
  
  return newSparseMatrices;
}

/////////////////////////////////////////////////////////////////
// Relax()
//
// Computes the consistency transformation for a single sequence
// z, and adds the transformed matrix to "posterior".
/////////////////////////////////////////////////////////////////

void AMAP_adapter::Relax (SparseMatrix *matXZ, SparseMatrix *matZY, VF &posterior){

  assert (matXZ);
  assert (matZY);

  int lengthX = matXZ->GetSeq1Length();
  int lengthY = matZY->GetSeq2Length();
  assert (matXZ->GetSeq2Length() == matZY->GetSeq1Length());

  // for every x[i]
  for (int i = 1; i <= lengthX; i++){
    SafeVector<PIF>::iterator XZptr = matXZ->GetRowPtr(i);
    SafeVector<PIF>::iterator XZend = XZptr + matXZ->GetRowSize(i);

    VF::iterator base = posterior.begin() + i * (lengthY + 1);

    // iterate through all x[i]-z[k]
    while (XZptr != XZend){
      SafeVector<PIF>::iterator ZYptr = matZY->GetRowPtr(XZptr->first);
      SafeVector<PIF>::iterator ZYend = ZYptr + matZY->GetRowSize(XZptr->first);
      const float XZval = XZptr->second;

      // iterate through all z[k]-y[j]
      while (ZYptr != ZYend){
        base[ZYptr->first] += XZval * ZYptr->second;
        ZYptr++;
      }
      XZptr++;
    }
  }

}


