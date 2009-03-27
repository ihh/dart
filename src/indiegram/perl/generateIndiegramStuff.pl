#!/usr/bin/perl -w


# I wrote this script to create some arrays, etc. which I needed for state typing and all that in my indiegram/tripletscfg.cc program.

# Used in SCFG_state_type_enum::State_type {}.
# List of state types, per the hash (xl yl zl zr yr xr).  (See tripletscfg.h for details.)
$stateTypeList = "";

# Used in SCFG_state_type_enum::state_type_string (State_type t) to get the string
# representation of an integer state type.
$stateTypeString = "";

# Used in SCFG_state_typing::_emit_size[] to get the emit size associated
# with a given state type.
$emitSize = "";

# Used in SCFG_state_typing::emit_xl_mul[], etc.
# Basically these are used to create a map between emissions (emit strings) and integers
# for array indexing.
$emitxlmul = "";
$emitxrmul = "";
$emitylmul = "";
$emityrmul = "";
$emitzlmul = "";
$emitzrmul = "";

# Now loop over all possible emissions.
# It's important, of course, to loop in the order 
# corresponding to the desired hash (xl yl zl zr yr xr).
for ($xl = 0; $xl < 2; ++$xl) {
  for ($yl = 0; $yl < 2; ++$yl) {
    for ($zl = 0; $zl < 2; ++$zl) {
      for ($zr = 0; $zr < 2; ++$zr) {
	for ($yr = 0; $yr < 2; ++$yr) {
	  for ($xr = 0; $xr < 2; ++$xr) {

	    # stateTypeList and stateTypeString
	    $stateName = "Emit";
	    if ($xl || $xr) {$stateName .= "X";}
	    if ($xl) {$stateName  .= "L";}
	    if ($xr) {$stateName .= "R";}
	    if ($yl || $yr) {$stateName .= "Y";}
	    if ($yl) {$stateName  .= "L";}
	    if ($yr) {$stateName .= "R";}
	    if ($zl || $zr) {$stateName .= "Z";}
	    if ($zl) {$stateName  .= "L";}
	    if ($zr) {$stateName .= "R";}
	    $state = $xl * 32 + $yl * 16 + $zl * 8 + $zr * 4 + $yr * 2 + $xr * 1;

	    $stateTypeList .= "$stateName = $state, ";
	    $stateTypeString .= "case $stateName: return \"$stateName\"; break;\n";

	    # emitSize
	    $numEmissions = $xl + $yl + $zl + $xr + $yr + $zr;
	    if ($numEmissions == 0) { $emitSize .= "_emit_N, "; } # because _emit_N = 0
	    else { $emitSize .= "_emit_$numEmissions, "; }

	    # emitxlmul, etc.
	    # remember that we're hashing as (xl yl zl zr yr xr)
	    # xr
	    if ($xr) {
	      $tmp = (0);
	    } else { $tmp = "N"; }
	    $emitxrmul .= "_emit_$tmp, ";

	    # yr
	    if ($yr) {
	      $tmp = ($xr);
	    } else { $tmp = "N"; }
	    $emityrmul .= "_emit_$tmp, ";

	    # zr
	    if ($zr) {
	      $tmp = ($xr+$yr);
	    } else { $tmp = "N"; }
	    $emitzrmul .= "_emit_$tmp, ";

	    # zl
	    if ($zl) {
	      $tmp = ($xr+$yr+$zr);
	    } else { $tmp = "N"; }
	    $emitzlmul .= "_emit_$tmp, ";

	    # yl
	    if ($yl) {
	      $tmp = ($xr+$yr+$zr+$zl);
	    } else { $tmp = "N"; }
	    $emitylmul .= "_emit_$tmp, ";

	    # xl
	    if ($xl) {
	      $tmp = ($xr+$yr+$zr+$zl+$yl);
	    } else { $tmp = "N"; }
	    $emitxlmul .= "_emit_$tmp, ";
	    
	  }
	}
	# add newlines for tidyness
	$stateTypeList .= "\n";
	$emitSize .= "\n";
	$emitxlmul .= "\n";
	$emitxrmul .= "\n";
	$emitylmul .= "\n";
	$emityrmul .= "\n";
	$emitzlmul .= "\n";
	$emitzrmul .= "\n";
      }
    }
  }
}

# now show the results!
print "SCFG_state_type_enum::State_type {}:\n";
print "$stateTypeList\n";

print "SCFG_state_type_enum::state_type_string (State_type t):\n";
print "$stateTypeString\n";

print "int SCFG_state_typing::_emit_size[64] = \n";
print "{ $emitSize\n";

print "int SCFG_state_typing::_emit_xl_mul[64] =\n";
print "{ $emitxlmul\n";

print "int SCFG_state_typing::_emit_xr_mul[64] =\n";
print "{ $emitxrmul\n";

print "int SCFG_state_typing::_emit_yl_mul[64] =\n";
print "{ $emitylmul\n";

print "int SCFG_state_typing::_emit_yr_mul[64] =\n";
print "{ $emityrmul\n";

print "int SCFG_state_typing::_emit_zl_mul[64] =\n";
print "{ $emitzlmul\n";

print "int SCFG_state_typing::_emit_zr_mul[64] =\n";
print "{ $emitzrmul\n";

