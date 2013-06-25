%fragment("Danu_NumPy_Fragments", "header",
          fragment="Danu_NumPy_Macros")


/* Given a PyObject return pointer to contiguous data */
%fragment("Danu_NumPy_Macros", "header")
{
%#define obj_is_float(a)         (PyFloat_Check(a))
%#define obj_is_int(a)           (PyInt_Check(a))
%#define obj_is_string(a)        (PyString_Check(a))

}
