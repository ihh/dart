## This is a concatenation of the .m4 autoconf macro definition files for Guile and pkg-config.
## I bundled them together as one file "acinclude.m4" because otherwise I couldn't get aclocal to load them individually...
## Ian Holmes, 2/26/2010

## Autoconf macros for working with Guile.
##
##   Copyright (C) 1998,2001, 2006 Free Software Foundation, Inc.
##
## This library is free software; you can redistribute it and/or
## modify it under the terms of the GNU Lesser General Public License
## as published by the Free Software Foundation; either version 3 of
## the License, or (at your option) any later version.
## 
## This library is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## Lesser General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public
## License along with this library; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301 USA

# serial 9

## Index
## -----
##
## GUILE_PROGS -- set paths to Guile interpreter, config and tool programs
## GUILE_FLAGS -- set flags for compiling and linking with Guile
## GUILE_SITE_DIR -- find path to Guile "site" directory
## GUILE_CHECK -- evaluate Guile Scheme code and capture the return value
## GUILE_MODULE_CHECK -- check feature of a Guile Scheme module
## GUILE_MODULE_AVAILABLE -- check availability of a Guile Scheme module
## GUILE_MODULE_REQUIRED -- fail if a Guile Scheme module is unavailable
## GUILE_MODULE_EXPORTS -- check if a module exports a variable
## GUILE_MODULE_REQUIRED_EXPORT -- fail if a module doesn't export a variable

## Code
## ----

## NOTE: Comments preceding an AC_DEFUN (starting from "Usage:") are massaged
## into doc/ref/autoconf-macros.texi (see Makefile.am in that directory).

# GUILE_PROGS -- set paths to Guile interpreter, config and tool programs
#
# Usage: GUILE_PROGS
#
# This macro looks for programs @code{guile}, @code{guile-config} and
# @code{guile-tools}, and sets variables @var{GUILE}, @var{GUILE_CONFIG} and
# @var{GUILE_TOOLS}, to their paths, respectively.  If either of the first two
# is not found, signal error.
#
# The variables are marked for substitution, as by @code{AC_SUBST}.
#
AC_DEFUN([GUILE_PROGS],
 [AC_PATH_PROG(GUILE,guile)
  if test "$GUILE" = "" ; then
# changed AC_MSG_ERROR to AC_WARN to allow compilation without guile - IH, 2/24/2010
      AC_WARN([guile not found - Scheme support will not be enabled])
      AC_SUBST(GUILE_INCLUDED,[0])
      AC_SUBST(GUILE_CDEFS)
  else
    AC_SUBST(GUILE)
    AC_PATH_PROG(GUILE_CONFIG,guile-config)
    if test "$GUILE_CONFIG" = "" ; then
# changed AC_MSG_ERROR to AC_WARN to allow compilation without guile - IH, 2/24/2010
      AC_WARN([guile-config not found - Scheme support will not be enabled])
      AC_SUBST(GUILE_INCLUDED,[0])
      AC_SUBST(GUILE_CDEFS)
    else
      AC_SUBST(GUILE_CONFIG)
      AC_PATH_PROG(GUILE_TOOLS,guile-tools)
      AC_SUBST(GUILE_TOOLS)
      AC_SUBST(GUILE_INCLUDED,[1])
      AC_SUBST(GUILE_CDEFS,[-DGUILE_INCLUDED])
    fi
  fi
 ])

# GUILE_FLAGS -- set flags for compiling and linking with Guile
#
# Usage: GUILE_FLAGS
#
# This macro runs the @code{guile-config} script, installed with Guile, to
# find out where Guile's header files and libraries are installed.  It sets
# two variables, @var{GUILE_CFLAGS} and @var{GUILE_LDFLAGS}.
#
# @var{GUILE_CFLAGS}: flags to pass to a C or C++ compiler to build code that
# uses Guile header files.  This is almost always just a @code{-I} flag.
#
# @var{GUILE_LDFLAGS}: flags to pass to the linker to link a program against
# Guile.  This includes @code{-lguile} for the Guile library itself, any
# libraries that Guile itself requires (like -lqthreads), and so on.  It may
# also include a @code{-L} flag to tell the compiler where to find the
# libraries.
#
# The variables are marked for substitution, as by @code{AC_SUBST}.
#
AC_DEFUN([GUILE_FLAGS],
 [AC_REQUIRE([GUILE_PROGS])dnl
  AC_MSG_CHECKING([libguile compile flags])
  GUILE_CFLAGS="`$GUILE_CONFIG compile`"
  AC_MSG_RESULT([$GUILE_CFLAGS])
  AC_MSG_CHECKING([libguile link flags])
  GUILE_LDFLAGS="`$GUILE_CONFIG link`"
  AC_MSG_RESULT([$GUILE_LDFLAGS])
  AC_SUBST(GUILE_CFLAGS)
  AC_SUBST(GUILE_LDFLAGS)
 ])

# GUILE_SITE_DIR -- find path to Guile "site" directory
#
# Usage: GUILE_SITE_DIR
#
# This looks for Guile's "site" directory, usually something like
# PREFIX/share/guile/site, and sets var @var{GUILE_SITE} to the path.
# Note that the var name is different from the macro name.
#
# The variable is marked for substitution, as by @code{AC_SUBST}.
#
AC_DEFUN([GUILE_SITE_DIR],
 [AC_REQUIRE([GUILE_PROGS])dnl
  AC_MSG_CHECKING(for Guile site directory)
  GUILE_SITE=`[$GUILE_CONFIG] info sitedir`
  AC_MSG_RESULT($GUILE_SITE)
  AC_SUBST(GUILE_SITE)
 ])

# GUILE_CHECK -- evaluate Guile Scheme code and capture the return value
#
# Usage: GUILE_CHECK_RETVAL(var,check)
#
# @var{var} is a shell variable name to be set to the return value.
# @var{check} is a Guile Scheme expression, evaluated with "$GUILE -c", and
#    returning either 0 or non-#f to indicate the check passed.
#    Non-0 number or #f indicates failure.
#    Avoid using the character "#" since that confuses autoconf.
#
AC_DEFUN([GUILE_CHECK],
 [AC_REQUIRE([GUILE_PROGS])
  $GUILE -c "$2" > /dev/null 2>&1
  $1=$?
 ])

# GUILE_MODULE_CHECK -- check feature of a Guile Scheme module
#
# Usage: GUILE_MODULE_CHECK(var,module,featuretest,description)
#
# @var{var} is a shell variable name to be set to "yes" or "no".
# @var{module} is a list of symbols, like: (ice-9 common-list).
# @var{featuretest} is an expression acceptable to GUILE_CHECK, q.v.
# @var{description} is a present-tense verb phrase (passed to AC_MSG_CHECKING).
#
AC_DEFUN([GUILE_MODULE_CHECK],
         [AC_MSG_CHECKING([if $2 $4])
	  GUILE_CHECK($1,(use-modules $2) (exit ((lambda () $3))))
	  if test "$$1" = "0" ; then $1=yes ; else $1=no ; fi
          AC_MSG_RESULT($$1)
         ])

# GUILE_MODULE_AVAILABLE -- check availability of a Guile Scheme module
#
# Usage: GUILE_MODULE_AVAILABLE(var,module)
#
# @var{var} is a shell variable name to be set to "yes" or "no".
# @var{module} is a list of symbols, like: (ice-9 common-list).
#
AC_DEFUN([GUILE_MODULE_AVAILABLE],
         [GUILE_MODULE_CHECK($1,$2,0,is available)
         ])

# GUILE_MODULE_REQUIRED -- fail if a Guile Scheme module is unavailable
#
# Usage: GUILE_MODULE_REQUIRED(symlist)
#
# @var{symlist} is a list of symbols, WITHOUT surrounding parens,
# like: ice-9 common-list.
#
AC_DEFUN([GUILE_MODULE_REQUIRED],
         [GUILE_MODULE_AVAILABLE(ac_guile_module_required, ($1))
          if test "$ac_guile_module_required" = "no" ; then
              AC_MSG_ERROR([required guile module not found: ($1)])
          fi
         ])

# GUILE_MODULE_EXPORTS -- check if a module exports a variable
#
# Usage: GUILE_MODULE_EXPORTS(var,module,modvar)
#
# @var{var} is a shell variable to be set to "yes" or "no".
# @var{module} is a list of symbols, like: (ice-9 common-list).
# @var{modvar} is the Guile Scheme variable to check.
#
AC_DEFUN([GUILE_MODULE_EXPORTS],
 [GUILE_MODULE_CHECK($1,$2,$3,exports `$3')
 ])

# GUILE_MODULE_REQUIRED_EXPORT -- fail if a module doesn't export a variable
#
# Usage: GUILE_MODULE_REQUIRED_EXPORT(module,modvar)
#
# @var{module} is a list of symbols, like: (ice-9 common-list).
# @var{modvar} is the Guile Scheme variable to check.
#
AC_DEFUN([GUILE_MODULE_REQUIRED_EXPORT],
 [GUILE_MODULE_EXPORTS(guile_module_required_export,$1,$2)
  if test "$guile_module_required_export" = "no" ; then
      AC_MSG_ERROR([module $1 does not export $2; required])
  fi
 ])

## guile.m4 ends here


# pkg.m4 - Macros to locate and utilise pkg-config.            -*- Autoconf -*-
# 
# Copyright © 2004 Scott James Remnant <scott@netsplit.com>.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
#
# As a special exception to the GNU General Public License, if you
# distribute this file as part of a program that contains a
# configuration script generated by Autoconf, you may include it under
# the same distribution terms that you use for the rest of that program.

# PKG_PROG_PKG_CONFIG([MIN-VERSION])
# ----------------------------------
AC_DEFUN([PKG_PROG_PKG_CONFIG],
[m4_pattern_forbid([^_?PKG_[A-Z_]+$])
m4_pattern_allow([^PKG_CONFIG(_PATH)?$])
AC_ARG_VAR([PKG_CONFIG], [path to pkg-config utility])dnl
if test "x$ac_cv_env_PKG_CONFIG_set" != "xset"; then
	AC_PATH_TOOL([PKG_CONFIG], [pkg-config])
fi
if test -n "$PKG_CONFIG"; then
	_pkg_min_version=m4_default([$1], [0.9.0])
	AC_MSG_CHECKING([pkg-config is at least version $_pkg_min_version])
	if $PKG_CONFIG --atleast-pkgconfig-version $_pkg_min_version; then
		AC_MSG_RESULT([yes])
	else
		AC_MSG_RESULT([no])
		PKG_CONFIG=""
	fi
		
fi[]dnl
])# PKG_PROG_PKG_CONFIG

# PKG_CHECK_EXISTS(MODULES, [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
#
# Check to see whether a particular set of modules exists.  Similar
# to PKG_CHECK_MODULES(), but does not set variables or print errors.
#
#
# Similar to PKG_CHECK_MODULES, make sure that the first instance of
# this or PKG_CHECK_MODULES is called, or make sure to call
# PKG_CHECK_EXISTS manually
# --------------------------------------------------------------
AC_DEFUN([PKG_CHECK_EXISTS],
[AC_REQUIRE([PKG_PROG_PKG_CONFIG])dnl
if test -n "$PKG_CONFIG" && \
    AC_RUN_LOG([$PKG_CONFIG --exists --print-errors "$1"]); then
  m4_ifval([$2], [$2], [:])
m4_ifvaln([$3], [else
  $3])dnl
fi])


# _PKG_CONFIG([VARIABLE], [COMMAND], [MODULES])
# ---------------------------------------------
m4_define([_PKG_CONFIG],
[if test -n "$$1"; then
    pkg_cv_[]$1="$$1"
 elif test -n "$PKG_CONFIG"; then
    PKG_CHECK_EXISTS([$3],
                     [pkg_cv_[]$1=`$PKG_CONFIG --[]$2 "$3" 2>/dev/null`],
		     [pkg_failed=yes])
 else
    pkg_failed=untried
fi[]dnl
])# _PKG_CONFIG

# _PKG_SHORT_ERRORS_SUPPORTED
# -----------------------------
AC_DEFUN([_PKG_SHORT_ERRORS_SUPPORTED],
[AC_REQUIRE([PKG_PROG_PKG_CONFIG])
if $PKG_CONFIG --atleast-pkgconfig-version 0.20; then
        _pkg_short_errors_supported=yes
else
        _pkg_short_errors_supported=no
fi[]dnl
])# _PKG_SHORT_ERRORS_SUPPORTED


# PKG_CHECK_MODULES(VARIABLE-PREFIX, MODULES, [ACTION-IF-FOUND],
# [ACTION-IF-NOT-FOUND])
#
#
# Note that if there is a possibility the first call to
# PKG_CHECK_MODULES might not happen, you should be sure to include an
# explicit call to PKG_PROG_PKG_CONFIG in your configure.ac
#
#
# --------------------------------------------------------------
AC_DEFUN([PKG_CHECK_MODULES],
[AC_REQUIRE([PKG_PROG_PKG_CONFIG])dnl
AC_ARG_VAR([$1][_CFLAGS], [C compiler flags for $1, overriding pkg-config])dnl
AC_ARG_VAR([$1][_LIBS], [linker flags for $1, overriding pkg-config])dnl

pkg_failed=no
AC_MSG_CHECKING([for $1])

_PKG_CONFIG([$1][_CFLAGS], [cflags], [$2])
_PKG_CONFIG([$1][_LIBS], [libs], [$2])

m4_define([_PKG_TEXT], [Alternatively, you may set the environment variables $1[]_CFLAGS
and $1[]_LIBS to avoid the need to call pkg-config.
See the pkg-config man page for more details.])

if test $pkg_failed = yes; then
        _PKG_SHORT_ERRORS_SUPPORTED
        if test $_pkg_short_errors_supported = yes; then
	        $1[]_PKG_ERRORS=`$PKG_CONFIG --short-errors --print-errors "$2" 2>&1`
        else 
	        $1[]_PKG_ERRORS=`$PKG_CONFIG --print-errors "$2" 2>&1`
        fi
	# Put the nasty error message in config.log where it belongs
	echo "$$1[]_PKG_ERRORS" >&AS_MESSAGE_LOG_FD

	ifelse([$4], , [AC_MSG_ERROR(dnl
[Package requirements ($2) were not met:

$$1_PKG_ERRORS

Consider adjusting the PKG_CONFIG_PATH environment variable if you
installed software in a non-standard prefix.

_PKG_TEXT
])],
		[AC_MSG_RESULT([no])
                $4])
elif test $pkg_failed = untried; then
	ifelse([$4], , [AC_MSG_FAILURE(dnl
[The pkg-config script could not be found or is too old.  Make sure it
is in your PATH or set the PKG_CONFIG environment variable to the full
path to pkg-config.

_PKG_TEXT

To get pkg-config, see <http://pkg-config.freedesktop.org/>.])],
		[$4])
else
	$1[]_CFLAGS=$pkg_cv_[]$1[]_CFLAGS
	$1[]_LIBS=$pkg_cv_[]$1[]_LIBS
        AC_MSG_RESULT([yes])
	ifelse([$3], , :, [$3])
fi[]dnl
])# PKG_CHECK_MODULES
# pkg.m4 ends here
