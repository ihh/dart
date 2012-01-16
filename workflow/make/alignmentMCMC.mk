
total_samples:=2000

%.handalign.stk: %.fa
	handalign $< -ts $(total_samples) -af $*.trace.stk -rs

# BAli-Phy will create its own directory for the files, labeled by run # and partition
# BAli-Phy is assumed in the user's path
%.baliphy: %.fa %.tree.nh
	bali-phy  $<  --tree $*.tree.nh -i $(total_samples)

.SECONDARY: