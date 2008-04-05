/* Need a small wrapper for the fortran main program, since fortran
   code should not be used as the main program, see information
   elsewhere */
extern "C" int mcotac1d();
int main()
{
    return mcotac1d();
}
