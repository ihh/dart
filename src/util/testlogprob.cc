#include "util/logprob.h"
#include "util/vector_output.h"

double tol = .01;
double InfProb = Score2Prob(InfinityScore);
bool approx (double a, double b)
{
  if (b >= InfProb && a >= InfProb) return TRUE;
  if (b <= -InfProb && a <= -InfProb) return TRUE;
  return abs(a-b)/(b==0?1:abs(b)) < tol;
}

// symbols with meaning: LogType, typestr, LVal1, LVal2

#define PASSOP(OP) cerr << "(" << typestr << ": " << LVal1 << " " << #OP << " " << LVal2 << " == " << (LVal1 OP LVal2) << " +/- " << (tol*100) << "%)\n";
#define PASSREL(REL) cerr << "(" << typestr << ": (" << LVal1 << " " << #REL << " " << LVal2 << ") is " << (LVal1 REL LVal2 ? "true" : "false") << ")\n";
#define FAILOP(OP,RVal,Target) THROWEXPR("Failed " << typestr << " test: (" << LVal1 << " " << #OP << " " << LVal2 << ") = " << RVal << ", should be " << Target << " +/- " << (tol*100) << "%)\n")
#define FAILREL(OP,RVal,Target) THROWEXPR("Failed " << typestr << " test: (" << LVal1 << " " << #OP << " " << LVal2 << ") = " << RVal << ", should be " << Target << "\n")

#define TestUnaryLL(OPEQ) { LogType x = LVal1; LogType y = LVal2; x OPEQ y; Prob p = LVal1; p OPEQ LVal2; if (!approx((double)x,p)) FAILOP(OPEQ,x,p); }
#define TestUnaryLP(OPEQ) { LogType x = LVal1; x OPEQ LVal2; Prob p = LVal1; p OPEQ LVal2; if (!approx((double)x,p)) FAILOP(OPEQ,x,p); }
#define TestUnary(OPEQ)   TestUnaryLL(OPEQ); TestUnaryLP(OPEQ);

#define TestBinaryLL(OP)  { LogType x = LVal1; LogType y = LVal2; LogType rval = x OP y; double p = LVal1 OP LVal2; if (!approx((double)rval,p)) FAILOP(OP,rval,p); }
#define TestBinaryLP(OP)  { LogType x = LVal1; LogType rval = x OP LVal2; double p = LVal1 OP LVal2; if (!approx((double)rval,p)) FAILOP(OP,rval,p); }
#define TestBinaryPL(OP)  { LogType x = LVal2; LogType rval = LVal1 OP x; double p = LVal1 OP LVal2; if (!approx((double)rval,p)) FAILOP(OP,rval,p); }

#define TestBinary(OP) \
TestBinaryLL(OP); \
TestBinaryLP(OP); \
TestBinaryPL(OP);

#define TestOp(OP) \
TestUnary(OP ## =); \
TestBinary(OP);

template<class LogType>
void testAllOp (const LogType& x, const char* typestr, double LVal1, double LVal2)
{
  TestOp(+);
  TestOp(-);
  TestOp(*);
  if (!approx(LVal1,0) || !approx(LVal2,0))  // avoid 0/0
    { TestOp(/); }  // enclose the multiline macro in brackets so it doesn't fall off the if clause
}

template<class LogType>
void testAllOpFB (const LogType& x, const char* typestr, double LVal1, double LVal2)
{
  testAllOp<LogType> (x, typestr, LVal1, LVal2);
  testAllOp<LogType> (x, typestr, LVal2, LVal1);
}

template<class LogType>
void testAllOpPM (const LogType& x, const char* typestr, double LVal1, double LVal2)
{
  // avoid -0 errors
  const bool zero1 = approx(LVal1,0);
  const bool zero2 = approx(LVal2,0);
  testAllOpFB<LogType> (x, typestr, LVal1, LVal2);
  if (!zero2)
    testAllOpFB<LogType> (x, typestr, LVal1, -LVal2);
  if (!zero1)
    testAllOpFB<LogType> (x, typestr, -LVal1, LVal2);
  if (!zero1 && !zero2)
    testAllOpFB<LogType> (x, typestr, -LVal1, -LVal2);
}

#define TestRelationLL(REL)  { LogType x = LVal1; LogType y = LVal2; bool rval = x REL y; bool p = LVal1 REL LVal2; if (rval ^ p) FAILREL(OP,rval,p); }
#define TestRelationLP(REL)  { LogType x = LVal1; bool rval = x REL LVal2; bool p = LVal1 REL LVal2; if (rval ^ p) FAILREL(OP,rval,p); }
#define TestRelationPL(REL)  { LogType x = LVal2; bool rval = LVal1 REL x; bool p = LVal1 REL LVal2; if (rval ^ p) FAILREL(OP,rval,p); }

#define TestRelation(REL) \
TestRelationLL(REL); \
TestRelationLP(REL); \
TestRelationPL(REL);

template<class LogType>
void testAllRelation (const LogType& x, const char* typestr, double LVal1, double LVal2)
{
  TestRelation(<);
  TestRelation(<=);
  TestRelation(>);
  TestRelation(>=);
  TestRelation(==);
  TestRelation(!=);
}

template<class LogType>
void testAllRelationFB (const LogType& x, const char* typestr, double LVal1, double LVal2)
{
  testAllRelation (x, typestr, LVal1, LVal2);
  testAllRelation (x, typestr, LVal2, LVal1);
  testAllRelation (x, typestr, LVal1, LVal1);
  testAllRelation (x, typestr, LVal2, LVal2);
}

template<class LogType>
void testAllRelationPM (const LogType& x, const char* typestr, double LVal1, double LVal2)
{
  testAllRelation (x, typestr, LVal1, LVal2);
  testAllRelation (x, typestr, LVal1, -LVal2);
  testAllRelation (x, typestr, -LVal1, LVal2);
  testAllRelation (x, typestr, -LVal1, -LVal2);
}

#define TestType(LogType,LVal1,LVal2) testAllOp (LogType(), #LogType, LVal1, LVal2); testAllRelationFB (LogType(), #LogType, LVal1, LVal2); cerr << "(tested " << LVal1 << " {+,-,*,/,<,<=,>,>=,==,!=} " << LVal2 << " for " << #LogType << ")\n";

#define TestTypePM(LogType,LVal1,LVal2) testAllOpPM (LogType(), #LogType, LVal1, LVal2); testAllRelationPM (LogType(), #LogType, LVal1, LVal2); cerr << "(tested {+,-}" << LVal1 << " {+,-,*,/,<,<=,>,>=,==,!=} {+,-}" << LVal2 << " for " << #LogType << ")\n";

void test_unsigned_types (Prob a, Prob b)
{
  if (a < b) THROWEXPR ("now this just ain't gonna work: a (" << a << ") should be greater than b (" << b << ")");

  tol = .01;
  TestType (Score_prob, a, b);

  tol = .001;
  TestType (Loge_prob, a, b);
  TestType (Slow_loge_prob, a, b);
}

void test_signed_types (Prob a, Prob b)
{
  tol = .01;
  TestTypePM (Score_real, a, b);

  tol = .001;
  TestTypePM (Loge_real, a, b);
  TestTypePM (Slow_loge_real, a, b);
}

void test_types (Prob a, Prob b)
{
  test_unsigned_types (a, b);
  test_signed_types (a, b);
}

#define TestCast(LogType,Val)    { double p = Val; LogType x = p; if (!approx ((double) x, p)) THROWEXPR ("Failed cast conversion for " << #LogType << ": returned " << x << ", expected " << p << "\n"); }
#define TestCastPM(LogType,Val)  TestCast(LogType,Val); TestCast(LogType,-Val);

void test_infinity()
{
  // unsigned types
  TestCast (Score_prob, InfProb);
  TestCast (Loge_prob, InfProb);
  TestCast (Slow_loge_prob, InfProb);
  // signed types
  TestCastPM (Score_real, InfProb);
  TestCastPM (Loge_real, InfProb);
  TestCastPM (Slow_loge_real, InfProb);
  cerr << "(tested infinity conversions)\n";
}

int main (int argc, char** argv)
{
  // parse command line
  INIT_OPTS_LIST (opts, argc, argv, 0, "[options]",
		  "logspace probability tester");
  opts.parse_or_die();

  // main code
  try
    {
      // give the arithmetic a workout
      test_types (0.4, 0.3);
      test_types (20, 10);
      test_types (1, 0);

      // test infinity conversions
      test_infinity();
    }
  catch (const Dart_exception& e)
    {
      CLOGERR << e.what();
      cout << "not ok\n";
      exit(1);
    }

  cerr << "(all score details taken care of. the bedrock of modern civilisation is safe)\n";
  cout << "ok\n";
  return 0;
}
