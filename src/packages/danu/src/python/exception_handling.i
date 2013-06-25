
%exception {
    clear_exception();
    $action;
    if (check_exception()) {
        PyErr_SetString(PyExc_RuntimeError, check_exception());
        return NULL;
    }
}

