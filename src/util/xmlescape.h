/* From http://mdawson.net/misc/xmlescape.php */

#include <string>
#include <sstream>
using namespace std;

/**
 * Escape characters that will interfere with xml.
 *
 * @param sSrc The src string to escape.
 * @return sSrc encoded for insertion into xml.
 */
std::string encodeForXml( const std::string &sSrc )
{
    ostringstream sRet;

    for( string::const_iterator iter = sSrc.begin(); iter!=sSrc.end(); iter++ )
    {
         unsigned char c = (unsigned char)*iter;

         switch( c )
         {
             case '&': sRet << "&amp;"; break;
             case '<': sRet << "&lt;"; break;
             case '>': sRet << "&gt;"; break;
             case '"': sRet << "&quot;"; break;
             case '\'': sRet << "&apos;"; break;

             default:
              if ( c<32 || c>127 )
              {
                   sRet << "&#" << (unsigned int)c << ";";
              }
              else
              {
                   sRet << c;
              }
         }
    }

    return sRet.str();
}
