  // Bifurcation iterators

  // Bifurc_in_iterator
  class Bifurc_in_iterator : public iterator <forward_iterator_tag, Bifurc_in> {
  private:
    // private data
    // left subseq:  by_start[N][lsubseq_idx]
    // right subseq: rsubseq
    Subseq_pointer subseq_begin, rsubseq, rsubseq_end;
    Subseq_index_rev_pointer lsubseq_idx, lsubseq_idx_end;
    Bifurc_in b;

    // private methods
    inline const Subseq& lsubseq() const { return subseq_begin[*lsubseq_idx]; }

  public:

    // end method
    inline Bifurc_in_iterator& goto_end() {
      rsubseq = rsubseq_end;
      lsubseq_idx = lsubseq_idx_end;
      return *this;
    }

    // init methods
    inline void init (Subseq_pointer sb,
	       Subseq_pointer r, Subseq_pointer re,
	       Subseq_index_rev_pointer l, Subseq_index_rev_pointer le)
    {
      subseq_begin = sb;
      rsubseq = r;
      rsubseq_end = re;
      lsubseq_idx = l;
      lsubseq_idx_end = le;
    }

    inline void init_copy (const Bifurc_in_iterator& b)
    {
      init (b.subseq_begin, b.rsubseq, b.rsubseq_end, b.lsubseq_idx, b.lsubseq_idx_end);
    }

    inline void init_full (const Fold_envelope_base& env, int subseq_idx)
    {
      const Subseq& subseq (env.subseq[subseq_idx]);
      const int s = subseq.start;
      const int e = subseq.end();
      const vector<int>& bs (env.by_start[s]);
      init (env.subseq.begin(),
	    env.by_end_begin(e), env.subseq.begin() + subseq_idx + 1,
	    bs.rbegin() + (bs.size() - subseq.by_start_index - 1), bs.rend());
    }

    // current element accessor
    inline const Bifurc_in* operator->() const { return &(operator*()); }
    inline const Bifurc_in& operator*() const {
      ((Bifurc_in&)b).l = *lsubseq_idx;
      ((Bifurc_in&)b).r = rsubseq - subseq_begin;
      return b;
    }

    // Provide both STL-like "==" & end() methods and Java-like finished() method
    inline const bool operator== (const Bifurc_in_iterator& i) const { return rsubseq == i.rsubseq; }
    inline const bool operator!= (const Bifurc_in_iterator& i) const { return rsubseq != i.rsubseq; }
    inline bool finished() const { return rsubseq == rsubseq_end; }

    // advance method
    inline Bifurc_in_iterator& operator++() {
      if (finished())
	THROWEXPR ("Attempt to advance iterator past end of range");
      do {
	// advance right subseq iterator (i.e. move its startpoint left)
	while (rsubseq->start >= lsubseq().end())
	  if (++rsubseq == rsubseq_end)  // right iterator at end of range?
	    return goto_end();
	// retreat left subseq iterator (i.e. move its endpoint left)
	while (lsubseq().end() > rsubseq->start)
	  if (++lsubseq_idx == lsubseq_idx_end)  // left iterator at end of range?
	    return goto_end();
      }	while (rsubseq->start != lsubseq().end());
      return *this;
    }
  };

  // Bifurc_outl_iterator
  class Bifurc_outl_iterator : public iterator <forward_iterator_tag, Bifurc_out_l> {
  private:
    // private data
    Subseq_pointer subseq_begin, osubseq, osubseq_end, lsubseq, lsubseq_end;
    Bifurc_out_l b;

  public:

    // end method
    inline Bifurc_outl_iterator& goto_end() {
      osubseq = osubseq_end;
      lsubseq = lsubseq_end;
      return *this;
    }

    // init methods
    inline void init (Subseq_pointer sb,
	       Subseq_pointer o, Subseq_pointer oe,
	       Subseq_pointer l, Subseq_pointer le)
    {
      subseq_begin = sb;
      osubseq = o;
      osubseq_end = oe;
      lsubseq = l;
      lsubseq_end = le;
    }

    inline void init_copy (const Bifurc_outl_iterator& b)
    {
      init (b.subseq_begin, b.osubseq, b.osubseq_end, b.lsubseq, b.lsubseq_end);
    }

    inline void init_full (const Fold_envelope_base& env, int subseq_idx)
    {
      const Subseq& subseq (env.subseq[subseq_idx]);
      const int s = subseq.start;
      const int e = subseq.end();
      init (env.subseq.begin(),
	    env.subseq.begin() + subseq_idx, env.by_end_end(e),
	    env.by_end_begin(s), env.by_end_end(s));
    }

    // current element accessor
    inline const Bifurc_out_l* operator->() const { return &(operator*()); }
    inline const Bifurc_out_l& operator*() const {
      ((Bifurc_out_l&)b).out = osubseq - subseq_begin;
      ((Bifurc_out_l&)b).l   = lsubseq - subseq_begin;
      return b;
    }

    // Provide both STL-like "==" & end() methods and Java-like finished()  method
    inline const bool operator== (const Bifurc_outl_iterator& i) const { return osubseq == i.osubseq; }
    inline const bool operator!= (const Bifurc_outl_iterator& i) const { return osubseq != i.osubseq; }
    inline bool finished() const { return osubseq == osubseq_end; }

    // advance method
    inline Bifurc_outl_iterator& operator++() {
      if (finished())
	THROWEXPR ("Attempt to advance iterator past end of range");
      do {
	// advance outside subseq iterator (i.e. move its startpoint left)
	while (osubseq->start >= lsubseq->start)
	  if (++osubseq == osubseq_end)  // outside iterator at end of range?
	    return goto_end();
	// advance left subseq iterator (i.e. move its startpoint left)
	while (lsubseq->start > osubseq->start)
	  if (++lsubseq == lsubseq_end)  // left iterator at end of range?
	    return goto_end();
	} while (osubseq->start != lsubseq->start);
      return *this;
    }
  };

  // Bifurc_outr_iterator
  class Bifurc_outr_iterator : public iterator <forward_iterator_tag, Bifurc_out_r> {
  private:
    // private data
    Subseq_pointer subseq_begin;
    Subseq_index_pointer osubseq_idx, osubseq_idx_end;
    Subseq_index_pointer rsubseq_idx, rsubseq_idx_end;
    Bifurc_out_r b;

    // private methods
    inline const Subseq& osubseq() const { return subseq_begin[*osubseq_idx]; }
    inline const Subseq& rsubseq() const { return subseq_begin[*rsubseq_idx]; }

  public:

    // end method
    inline Bifurc_outr_iterator& goto_end() {
      osubseq_idx = osubseq_idx_end;
      rsubseq_idx = rsubseq_idx_end;
      return *this;
    }

    // init methods
    inline void init (Subseq_pointer sb,
	       Subseq_index_pointer o, Subseq_index_pointer oe,
	       Subseq_index_pointer r, Subseq_index_pointer re)
    {
      subseq_begin = sb;
      osubseq_idx = o;
      osubseq_idx_end = oe;
      rsubseq_idx = r;
      rsubseq_idx_end = re;
    }

    inline void init_copy (const Bifurc_outr_iterator& b)
    {
      init (b.subseq_begin, b.osubseq_idx, b.osubseq_idx_end, b.rsubseq_idx, b.rsubseq_idx_end);
    }

    inline void init_full (const Fold_envelope_base& env, int subseq_idx)
    {
      const Subseq& subseq (env.subseq[subseq_idx]);
      const int s = subseq.start;
      const int e = subseq.end();
      init (env.subseq.begin(),
	    env.by_start[s].begin() + subseq.by_start_index, env.by_start[s].end(),
	    env.by_start[e].begin(), env.by_start[e].end());
    }

    // current element accessor
    inline const Bifurc_out_r* operator->() const { return &(operator*()); }
    inline const Bifurc_out_r& operator*() const {
      ((Bifurc_out_r&)b).out = *osubseq_idx;
      ((Bifurc_out_r&)b).r   = *rsubseq_idx;
      return b;
    }

    // Provide both STL-like "==" & end() methods and Java-like finished()  method
    inline const bool operator== (const Bifurc_outr_iterator& i) const { return osubseq_idx == i.osubseq_idx; }
    inline const bool operator!= (const Bifurc_outr_iterator& i) const { return osubseq_idx != i.osubseq_idx; }
    inline bool finished() const { return osubseq_idx == osubseq_idx_end; }

    // advance method
    inline Bifurc_outr_iterator& operator++() {
      if (finished())
	THROWEXPR ("Attempt to advance iterator past end of range");
      do {
	// advance outside subseq iterator (i.e. move its endpoint right)
	while (osubseq().end() <= rsubseq().end())
	  if (++osubseq_idx == osubseq_idx_end)  // outside iterator at end of range?
	    return goto_end();
	// advance right subseq iterator (i.e. move its endpoint right)
	while (rsubseq().end() < osubseq().end())
	  if (++rsubseq_idx == rsubseq_idx_end)  // right iterator at end of range?
	    return goto_end();
      } while (osubseq().end() != rsubseq().end());
      return *this;
    }
  };

  // Bifurcation iterator pseudo-vectors
  template <class Bifurc_iterator>
  struct Bifurc_pseudovec
  {
    // typedefs
    typedef Bifurc_iterator iterator;
    typedef Bifurc_iterator const_iterator;
    typedef typename Bifurc_iterator::value_type value_type;
    // data
    Bifurc_iterator b, e;
    vector<value_type> my_values;
    bool my_values_valid;
    // initialiser
    inline void init (const Fold_envelope_base& env, int subseq_idx)
    { my_values_valid = false; b.init_full (env, subseq_idx); e.init_copy (b); e.goto_end(); }
    // begin & end methods
    inline const iterator& begin() const { return b; }
    inline const iterator& end() const { return e; }
    // slow methods for backward compatibility with vectors - avoid
    // (these reconstruct the original vectors & are very inefficient)
    inline size_t size() const
    {
      if (!my_values_valid)
	build_my_values();
      return my_values.size();
    }
    inline value_type operator[] (size_t n) const
    {
      if (!my_values_valid)
	build_my_values();
      return my_values[n];
    }
    inline void build_my_values() const
    {
      Bifurc_pseudovec<iterator>* mutable_this = (Bifurc_pseudovec<iterator>*) this;
      mutable_this->my_values.clear();
      for_iterator (iterator, i, begin(), end())
	mutable_this->my_values.push_back (*i);
      mutable_this->my_values_valid = true;
    }
  };

#ifdef DART_USE_BIFURCATION_PSEUDOVECTORS
  typedef Bifurc_pseudovec<Bifurc_in_iterator> Bifurc_in_pseudovec;
  typedef Bifurc_pseudovec<Bifurc_outl_iterator> Bifurc_outl_pseudovec;
  typedef Bifurc_pseudovec<Bifurc_outr_iterator> Bifurc_outr_pseudovec;
#else /* DART_USE_BIFURCATION_PSEUDOVECTORS */
  typedef vector<Bifurc_in> Bifurc_in_pseudovec;
  typedef vector<Bifurc_out_l> Bifurc_outl_pseudovec;
  typedef vector<Bifurc_out_r> Bifurc_outr_pseudovec;
#endif /* DART_USE_BIFURCATION_PSEUDOVECTORS */

