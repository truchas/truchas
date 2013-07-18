#
# STRING_TO_LIST(string list)
#
# Convert whitespaced string into a list variable
#
MACRO(STRING_TO_LIST str newlist)
  string(REGEX MATCHALL "([^\ ]+\ |[^\ ]+$)" ${newlist} "${str}")
ENDMACRO(STRING_TO_LIST string list)
