#include "scfg/foldenv.h"
#include "indiegram/tripletscfgdp.h"
#include "indiegram/postprobs.h"

#define DOTPLOT_REPORT_INTERVAL 10000  /* number of cells between logfile messages during DP */

Triplet_inside_outside_matrix::Triplet_inside_outside_matrix (const Triplet_SCFG& scfg,
							      const Named_profile& np_x, const Named_profile& np_y, const Named_profile& np_z,
							      const Fold_envelope& foldenv_x, const Fold_envelope& foldenv_y, const Fold_envelope& foldenv_z,
							      bool fill_now /* = true */)
  : inside (scfg, np_x, np_y, np_z, foldenv_x, foldenv_y, foldenv_z, fill_now),
    outside (inside, fill_now)
{ }

void Triplet_inside_outside_matrix::show (ostream& o) const
{
  Biosequence xseq, yseq, zseq;
  inside.scfg.alphabet().dsq2seq (inside.dsq_x, xseq);
  inside.scfg.alphabet().dsq2seq (inside.dsq_y, yseq);
  inside.scfg.alphabet().dsq2seq (inside.dsq_z, zseq);
  o << "Sequence X: " << xseq << "\n";
  o << "Sequence Y: " << yseq << "\n";
  o << "Sequence Z: " << zseq << "\n";
  o << "Envelope X:\n";
  inside.foldenv_x.dump(o);
  o << "Envelope Y:\n";
  inside.foldenv_y.dump(o);
  o << "Envelope Z:\n";
  inside.foldenv_z.dump(o);
  o << "Inside matrix:\n";
  inside.show(o);
  o << "Outside matrix:\n";
  outside.show(o);

}



Score Triplet_inside_outside_matrix::post_transition_sc (int src_state, int dest_state, int subseq_idx_x, int subseq_idx_y, int subseq_idx_z) const
{
  // make copies of some variables
  const Triplet_SCFG& scfg = inside.scfg;
  const Subseq& subseq_x = inside.foldenv_x.subseq[subseq_idx_x];
  const Subseq& subseq_y = inside.foldenv_y.subseq[subseq_idx_y];
  const Subseq& subseq_z = inside.foldenv_z.subseq[subseq_idx_z];

  // calculate P(inside) / P(final)
  Score inside_minus_final_sc = -inside.final_sc;
  if (dest_state != Grammar_state_enum::End)
    ScorePMulAcc (inside_minus_final_sc, inside.read_cell (dest_state, subseq_idx_x, subseq_idx_y, subseq_idx_z));
  else  // dest_state == End
    if (subseq_x.len > 0 || subseq_y.len > 0 || subseq_z.len > 0)
      return -InfinityScore;

  // calculate P(outside) * P(transition)
  Score incoming_sc = scfg.transition_scores.transition (src_state, dest_state);
  if (src_state != Grammar_state_enum::Start)
    ScorePMulAcc (incoming_sc, outside.read_cell (src_state, subseq_idx_x, subseq_idx_y, subseq_idx_z));
  else  // src_state == Start
    if (!inside.local_start_allowed (subseq_x, subseq_y, subseq_z))
      return -InfinityScore;

  // multiply the probabilities and return
  return ScorePMul (incoming_sc, inside_minus_final_sc);
}

Score Triplet_inside_outside_matrix::post_state_sc (int dest_state, int subseq_idx_x, int subseq_idx_y, int subseq_idx_z) const
{
  // get subseqs
  const Subseq& subseq_x = inside.foldenv_x.subseq[subseq_idx_x];
  const Subseq& subseq_y = inside.foldenv_y.subseq[subseq_idx_y];
  const Subseq& subseq_z = inside.foldenv_z.subseq[subseq_idx_z];

  // sum posterior probabilities of incoming transitions
  Score sc = post_transition_sc (Grammar_state_enum::Start, dest_state, subseq_idx_x, subseq_idx_y, subseq_idx_z);
  for_const_contents (vector<int>, inside.allowed_src_out_states (dest_state, subseq_x, subseq_y, subseq_z), s)
    ScorePSumAcc (sc, post_transition_sc (*s, dest_state, subseq_idx_x, subseq_idx_y, subseq_idx_z));
  return sc;
}




void Triplet_fold_dotplot::write_dotplot (const sstring& filename) const
{
  ofstream file (filename.c_str());
  if (!file)
    THROWEXPR ("Couldn't create dotplot file with name '" << filename << "'");
  // print horizontal sequence axis labels
  file << ".";
  for (int x = 0; x < (int) seq.size(); ++x)
    file << ' ' << seq[x];
  file << '\n';
  // print rows
  for (int y = 0; y < (int) seq.size(); ++y) {
    file << seq[y];
    for (int x = 0; x < (int) seq.size(); ++x)
      file << ' ' << (*this) (x,y);
    file << '\n';
  }
}

void Triplet_fold_dotplot::write_dotplot (const sstring& prefix, const sstring& seqname) const
{
  sstring filename;
  filename << prefix << '-' << seqname;
  write_dotplot (filename);
}




Triplet_fold_dotplot::Triplet_fold_dotplot (const Triplet_inside_outside_matrix& in_out, int seq_index)
  : array2d<Prob> (seq_index == 0 ? in_out.inside.np_x.seq.size() : (seq_index == 1 ? in_out.inside.np_y.seq.size() : in_out.inside.np_z.seq.size()),
		   seq_index == 0 ? in_out.inside.np_x.seq.size() : (seq_index == 1 ? in_out.inside.np_y.seq.size() : in_out.inside.np_z.seq.size()), 0.),
    seq (seq_index == 0 ? in_out.inside.np_x.seq : (seq_index == 1 ? in_out.inside.np_y.seq : in_out.inside.np_z.seq))
{
  const Fold_envelope& foldenv_x (in_out.inside.foldenv_x);
  const Fold_envelope& foldenv_y (in_out.inside.foldenv_y);
  const Fold_envelope& foldenv_z (in_out.inside.foldenv_z);

  const Triplet_SCFG& scfg (in_out.inside.scfg);

  const Fold_envelope& foldenv (seq_index == 0 ? foldenv_x : (seq_index == 1 ? foldenv_y : foldenv_z));

  // foldenv: other_foldenv_1, other_foldenv_2
  // x: y, z
  // y: x, z
  // z: x, y
  const Fold_envelope& other_foldenv_1 (seq_index == 0 ? foldenv_y : foldenv_x);
  const Fold_envelope& other_foldenv_2 (seq_index == 2 ? foldenv_y : foldenv_z);

  const State_type lr_emit_flags = (seq_index == 0 ? EmitXLR : (seq_index == 1 ? EmitYLR : EmitZLR));

  const int seqlen = seq.size();

  CTAG(4,DOTPLOT) << "Starting fold dotplot prob calculations for sequence "
		  << (seq_index == 0 ? in_out.inside.np_x.name : (seq_index == 1 ? in_out.inside.np_y.name : in_out.inside.np_z.name)) << ".\n";
  const bool logging = CTAGGING(1,DOTPLOT);

  const int total_cells = seqlen * (seqlen + 1) / 2;
  int cells_done = 0, last_report = 0;

  for (int start = 0; start < seqlen; ++start)
    for (int end = start + 1; end <= seqlen; ++end) {
      const int subseq_idx = foldenv.find_subseq_idx (start, end - start);
      if (subseq_idx >= 0)
	{
	  Prob& p = (*this) (end - 1, start);
	  for (int s = 0; s < scfg.num_states(); ++s)
	    if ((scfg.state_type[s] & lr_emit_flags) == lr_emit_flags)
	      for (int other_subseq_idx_1 = 0; other_subseq_idx_1 < other_foldenv_1.subseqs(); ++other_subseq_idx_1) {
		for (int other_subseq_idx_2 = 0; other_subseq_idx_2 < other_foldenv_2.subseqs(); ++other_subseq_idx_2) {
		  int subseq_idx_x, subseq_idx_y, subseq_idx_z;
		  if (seq_index == 0) {
		    subseq_idx_x = subseq_idx;
		    subseq_idx_y = other_subseq_idx_1;
		    subseq_idx_z = other_subseq_idx_2;
		  } else if (seq_index == 1) {
		    subseq_idx_x = other_subseq_idx_1;
		    subseq_idx_y = subseq_idx;
		    subseq_idx_z = other_subseq_idx_2;
		  } else {
		    subseq_idx_x = other_subseq_idx_1;
		    subseq_idx_y = other_subseq_idx_2;
		    subseq_idx_z = subseq_idx;
		  }
		  const Prob q = Score2Prob (in_out.post_state_sc (s, subseq_idx_x, subseq_idx_y, subseq_idx_z));
		  p += q;
		  if (logging)
		    CTAG(1,DOTPLOT) << "P (Emit" << (seq_index==0 ? "X" : (seq_index==1 ? "Y" : "Z")) << "LR,state=" << s << ",subseq_x=[" << foldenv_x.subseq[subseq_idx_x].start << ".." << foldenv_x.subseq[subseq_idx_x].end() << "],subseq_y=[" << foldenv_y.subseq[subseq_idx_y].start << ".." << foldenv_y.subseq[subseq_idx_y].end() << "],subseq_z=[" << foldenv_z.subseq[subseq_idx_z].start << ".." << foldenv_z.subseq[subseq_idx_z].end() << "]) = " << q << "\n";
		}
	      }

	  // print this entry to log
	  CTAG(2,DOTPLOT) << "P (" << scfg.state_type_string (lr_emit_flags) << ",start=" << start << ",end=" << end << ") = " << p << "\n";

	  ++cells_done;
	  if (cells_done - last_report >= DOTPLOT_REPORT_INTERVAL
	      || (cells_done == total_cells && cells_done > last_report))
	    {
	      last_report = cells_done;
	      CTAG(1,CFGDP) << "Triplet_fold_dotplot: finished " << cells_done << " basepairs (" << ((int)(1000.*(double)cells_done/(double)total_cells))/10. << "%)\n";
	    }
	}
    }
}



