#include "seq/suffix.h"
#include "util/logfile.h"
#include "util/vector_output.h"

Word_count::Word_count (const Sequence_database_index& index, int len) : word_len (len)
{
  CTAG(5,WORD_COUNT) << "Counting " << len << "-mer frequencies\n";
  vector<int> key (len);
  for (int seq = 0; seq < index.size(); ++seq)
    {
      const Digitized_biosequence& dsq = index.profile[seq]->dsq;
      const int seq_len = dsq.size();
      for (int pos = 0; pos < seq_len + 1 - len; ++pos)
	{
	  copy (dsq.begin() + pos, dsq.begin() + pos + len, key.begin());
	  ++(*this)[key];
	}
      if ((seq+1) % 500 == 0) CTAG(4,WORD_COUNT) << "Counted " << (seq+1) << " sequences\n";
    }
  if (index.size() % 500 != 0) CTAG(4,WORD_COUNT) << "Counted " << index.size() << " sequences\n";
  if (CTAGGING(3,WORD_COUNT))
    {
      CL << len << "-mer frequencies:\n";
      for_const_contents (Word_count, *this, wc)
	CL << "(" << wc->first << "): " << wc->second << "\n";
    }
}

int Word_count::prefix_count (const int* prefix_begin, const int* prefix_end) const
{
  if (prefix_end - prefix_begin > word_len) return 0;
  int n = 0;
  for_const_contents (Word_count, *this, wc)
    if (equal (prefix_begin, prefix_end, wc->first.begin()))
      n += wc->second;
  return n;
}

Suffix_tree::~Suffix_tree() { delete_nodes(); }

void Suffix_tree::clear()
{
  delete_nodes();
  node.clear();
}

void Suffix_tree::swap (Suffix_tree& tree)
{
  node.swap (tree.node);
}

void Suffix_tree::delete_nodes()
{
  for_contents (vector<Node*>, node, n) delete *n;
}

void Suffix_tree::add_suffix (const int* suff_begin, const int* suff_end)
{
  const int len = suff_end - suff_begin;
  int n = 0;  // current node
  int i = 0;  // current offset to suff_begin
  ++node[0]->count;  // increment count of root node
  while (i < len)
    {
      // look for exact or partial match
      bool exact = 0;
      bool partial = 0;
      for_contents (vector<int>, node[n]->child, iter)
	{
	  const vector<int>& key = node[*iter]->data;
	  int match = 0;
	  for (int j = 0; j < min ((int) key.size(), len-i); ++j)
	    if (key[j] == suff_begin[i+j]) ++match; else break;
	  if (match == (int) key.size())  // exact match
	    {
	      ++node[*iter]->count;
	      i += match;
	      n = *iter;
	      exact = 1;
	      // CL << "Matched up to " << i << "; going to node " << n << "\n";
	      break;
	    }
	  if (match > 0)  // partial match, so split the branch
	    {
	      const vector<int> key_base (key.begin(), key.begin() + match);
	      const vector<int> old_key (key.begin() + match, key.end());
	      const vector<int> new_key (suff_begin + i + match, suff_end);

	      const int split_node = node.size();
	      const int new_node = node.size() + 1;

	      const int old_node = *iter;
	      *iter = split_node;
	      
	      node.push_back (new Node());  // invalidates iter

	      node[split_node]->child.push_back (old_node);

	      node[split_node]->data = key_base;
	      node[old_node]->data = old_key;

	      node[split_node]->count = node[old_node]->count;
	      ++node[split_node]->count;
	      
	      if (i + match < len)  // don't add a new node if new suffix is a subsequence of old suffix
		{
		  node.push_back (new Node());
		  node[split_node]->child.push_back (new_node);
		  node[new_node]->data = new_key;
		  ++node[new_node]->count;
		}

	      // CL << "Adding node " << split_node << " before " << old_node;
	      // if (i + match < len) CL << " and new node " << new_node;
	      // CL << "; new tree:\n";
	      // dump (CL, 0);

	      partial = 1;
	      break;
	    }
	}
      if (exact) continue;
      if (partial) break;

      // make new node
      const int match_len = allow_long_branches ? (suff_end - suff_begin - i) : 1;
      const vector<int> new_key (suff_begin + i, suff_begin + i + match_len);
      node.push_back (new Node());
      const int new_node = node.size() - 1;
      node[n]->child.push_back (new_node);
      ++node[new_node]->count;
      node[new_node]->data = new_key;
      n = new_node;
      i += match_len;

      // CL << "Adding node " << new_node << "; new tree:\n";
      // dump (CL, 0);
    }
}

void Suffix_tree::add_suffices (const vector<int>& seq_data, int max_len)
{
  const int seq_len = seq_data.size();
  const int max_suff_len = max_len ? max_len : seq_len;
  for (int pos = 0; pos < (int) seq_data.size(); ++pos)
    add_suffix (&seq_data[pos], &seq_data[min (pos + max_suff_len, seq_len)]);
}

int Suffix_tree::find_node (const int* suff_begin, const int* suff_end)
{
  const int len = suff_end - suff_begin;
  int n = 0;  // current node
  int i = 0;  // current offset to suff_begin
  while (i < len)
    {
      // look for exact or partial match
      bool exact = 0;
      for_const_contents (vector<int>, node[n]->child, iter)
	{
	  const vector<int>& key = node[*iter]->data;
	  int match = 0;
	  for (int j = 0; j < min ((int) key.size(), len-i); ++j)
	    if (key[j] == suff_begin[i+j]) ++match; else break;
	  if (match == (int) key.size())  // exact match
	    {
	      i += match;
	      n = *iter;
	      exact = 1;
	      break;
	    }
	  if (match > 0) break;  // partial match, so stop
	}
      if (exact) continue;
      n = -1;
      break;
    }
  return n;
}

int Suffix_tree::suffix_count (const int* suff_begin, const int* suff_end)
{
  const int n = find_node (suff_begin, suff_end);
  return n < 0 ? 0 : node[n]->count;
}

Suffix_tree::Suffix_tree (bool allow_long_branches)
  : allow_long_branches (allow_long_branches),
    node (1, new Node())
{ }

void Suffix_tree::build (const Sequence_database_index& index, int max_len, int min_len, Word_count* word_count, int min_count)
{
  int chars = 0;
  for (int seq = 0; seq < index.size(); ++seq)
    {
      const Digitized_biosequence& dsq = index.profile[seq]->dsq;
      const int seq_len = dsq.size();
      const int max_suff_len = max_len ? max_len : seq_len;
      if (word_count)
	for (int pos = 0; pos < seq_len + 1 - min_len; ++pos)
	  {
	    const vector<int> seed (&dsq[pos], &dsq[pos + min_len]);
	    if ((*word_count)[seed] >= min_count)
	      add_suffix (&dsq[pos], &dsq[min (pos + max_suff_len, seq_len)]);
	  }
      else
	add_suffices (dsq, max_len);
      chars += dsq.size();
      if (CTAGGING(3,SUFFIX)) CL << "Counted sequence " << seq+1 << "/" << index.size() << " '" << index.profile[seq]->name << "'; " << node.size() << " nodes, " << chars << " characters\n";
      if ((seq+1) % 500 == 0) CTAG(4,SUFFIX) << "Counted " << seq+1 << " sequences\n";
    }
  if (index.size() % 500 != 0) CTAG(4,SUFFIX) << "Counted " << index.size() << " sequences\n";
}

Suffix_tree::Suffix_tree (const Sequence_database_index& index, int max_len, int min_len, int min_count, bool allow_long_branches)
  : allow_long_branches (allow_long_branches),
    node (1, new Node())
{
  if (max_len && max_len < min_len) return;
  if (min_count > 1 && min_len > 1)
    {
      Word_count word_count (index, min_len);
      CTAG(5,SUFFIX) << "Building suffix tree\n";
      build (index, max_len, min_len, &word_count, min_count);
    }
  else
    build (index, max_len, min_len, 0, 0);
}

void Suffix_tree::dump (ostream& o, const Alphabet* alph) const
{
  for (int n = 0; n < (int) node.size(); ++n)
    {
      Digitized_biosequence dsq = node[n]->data;
      if (alph)
	{
	  Biosequence seq;
	  alph->dsq2seq (dsq, seq);
	  o << n << " \"" << seq << "\"";
	}
      else
	o << "\"" << dsq << "\"";
      if (node[n]->child.size()) o << "->(" << node[n]->child << ")";
      o << ": " << node[n]->count << "\n";
    }
}

void Suffix_tree::set_null_metascores (const Sequence_database_index& index, int metasc_idx, int null_order, int alphabet_size)
{
  if (CTAGGING(3,SUFFIX)) CL << "Setting order-" << null_order << " null model metascores\n";
  for (int seq = 0; seq < index.size(); ++seq)
    {
      const Digitized_biosequence& dsq = index.profile[seq]->dsq;
      Metascore& meta = index.profile[seq]->meta_sc[metasc_idx];
      const int seq_len = dsq.size();
      for (int pos = 0; pos < seq_len; ++pos)
	{
	  int prev_node = 0;
	  for (int word_len = min(pos,null_order); word_len > 0; --word_len)  // get the maximum allowable context
	    {
	      int n = find_node (&dsq[pos - word_len], &dsq[pos]);
	      if (n >= 0) { prev_node = n; break; }
	    }
	  const int prev_count = node[prev_node]->count;
	  int this_node = -1;
	  for_const_contents (vector<int>, node[prev_node]->child, c)
	    if (node[*c]->data.size() == 1 && node[*c]->data[0] == dsq[pos])
	      { this_node = *c; break; }
	  const int this_count = this_node < 0 ? 0 : node[this_node]->count;
	  const Prob p = ((double) (this_count + 1)) / ((double) (prev_count + alphabet_size));  // NB Laplace pseudocounts
	  meta[pos] = Prob2Score (p);
	}
      if (CTAGGING(1,SUFFIX)) CL << "Null scores for '" << index.name[seq] << "': (" << meta << ")\n";
      if ((seq+1) % 500 == 0) CTAG(4,SUFFIX) << "Computed null scores for " << seq+1 << " sequences\n";
    }
  if (index.size() % 500 != 0) CTAG(4,SUFFIX) << "Computed null scores for " << index.size() << " sequences\n";
}
