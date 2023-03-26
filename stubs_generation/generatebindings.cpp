// Original author: Aziem Chawdhary
// Slightly adapted https://github.com/aziem/libclang-ocaml-bindings-generator

#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include "cmdparser.hpp"
#include <clang-c/Index.h>

using namespace std;

vector<string> blacklist = {"pari_fprintf",      "pari_fread_chars",
                            "pari_vfprintf",     "pari_fclose",
                            "gp_readvec_stream", "gp_read_stream",
                            "file_is_binary",    "pari_vfprintf",
                            "pari_vprintf",      "pari_vsprintf",
                            "out_vprintf",       "gvsprintf",
                            "os_signal",         "out_print0",
                            "out_print1",        "out_printf",
                            "out_putc",          "out_puts",
                            "out_term_color",    "pari_stack_base",
                            "pari_stack_init",   "pari_stack_alloc",
                            "pari_stack_new",    "pari_stack_delete",
                            "pari_stack_pushp",  "pari_stackcheck_init",
                            "pari_stack",        "PariOUT",
                            "switchin",          "switchout"};

struct EnumDecl {
  string name;
  vector<string> enumconstants;
};

class FuncDecl {
public:
  string name;
  string resulttype;
  vector<string> paramtypes;

  bool operator==(const FuncDecl &f) { return (this->name == f.name); }

  int operator<(const FuncDecl &f) { return (this->name < f.name); }

  int operator<(FuncDecl f) { return (this->name < f.name); }

  friend int operator<(const string s, const FuncDecl f);
};
int operator<(const FuncDecl f, string s) { return (f.name < s); }

int operator<(const string s, const FuncDecl f) { return (s < f.name); }

class StructDecl {
public:
  string name;
  vector<tuple<string, string>> fieldnames;
  bool is_typedef;

  bool operator==(const StructDecl &f) { return (this->name == f.name); }

  int operator<(const StructDecl &f) { return (this->name < f.name); }

  int operator<(StructDecl f) { return (this->name < f.name); }

  friend int operator<(const string s, const StructDecl f);
};
int operator<(const StructDecl f, string s) { return (f.name < s); }

int operator<(const string s, const StructDecl f) { return (s < f.name); }

struct TypedefDecl {
  string name;
  string type;
  string ocaml_type;
};

ostream &operator<<(ostream &stream, const CXString &str) {
  stream << clang_getCString(str);
  clang_disposeString(str);
  return stream;
}

string getFromCXString(const CXString &str) {
  string s(clang_getCString(str));
  clang_disposeString(str);
  return s;
}

string lowerCase(string s) {
  string lowername = s;
  transform(lowername.begin(), lowername.end(), lowername.begin(), ::tolower);
  return lowername;
}

string upperCase(string s) {
  string uppername = s;
  transform(uppername.begin(), uppername.end(), uppername.begin(), ::toupper);
  return uppername;
}

static inline void replaceAll(std::string &str, const std::string &from,
                              const std::string &to) {
  size_t start_pos = 0;
  while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
    str.replace(start_pos, from.length(), to);
    start_pos += to.length();
  }
}

CXTranslationUnit parseFile(CXIndex index, string filename, bool iscpp) {
  CXTranslationUnit unit;
  if (iscpp) {
    const char *const opt = "-xc++";
    unit = clang_parseTranslationUnit(index, filename.c_str(), &opt, 1, nullptr,
                                      0, CXTranslationUnit_None);
  } else {
    unit = clang_parseTranslationUnit(index, filename.c_str(), nullptr, 0,
                                      nullptr, 0, CXTranslationUnit_None);
  }
  if (unit == nullptr) {
    cerr << "Clang parsing failed." << endl;
    exit(-1);
  } else {
    return unit;
  }
}

string typeToStringOCaml(CXType ty, CXCursor c) {
  switch (ty.kind) {
  case CXType_Void:
    return "void";
  case CXType_UInt:
    return "size_t";
  case CXType_Int:
    return "int";
  case CXType_ULong:
    return "Unsigned.ULong.t";
  case CXType_Long:
    return "Signed.Long.t";
  case CXType_Bool:
    return "bool";
  case CXType_UChar:
    return "Unsigned.uchar";
  case CXType_Char_S:
  case CXType_Char_U:
    return "char";
  case CXType_ConstantArray: {
    string s = getFromCXString(clang_getTypeSpelling(ty));
    return "(" + typeToStringOCaml(clang_getArrayElementType(ty), c) +
           ") carray";
  }
  case CXType_IncompleteArray: {
    string s = getFromCXString(clang_getTypeSpelling(ty));
    return "(" + typeToStringOCaml(clang_getArrayElementType(ty), c) + ") ptr";
  }
  case CXType_Typedef: {
    string s = getFromCXString(clang_getTypeSpelling(ty));
    // Ugly hack: check if the typedef contains an int and if so, use
    // the OCaml equivalent.
    size_t found = s.find("int");
    if (found != std::string::npos) {
      return s;
    }
    // Otherwise return the typedef'd type.
    return getFromCXString(clang_getTypeSpelling(ty));
  }
  case CXType_FunctionProto: {
    // For function proto we need to do some extra work by traversing
    // the cursor to get the parameter types.
    string rettype = typeToStringOCaml(clang_getResultType(ty), c);
    vector<string> parmtypes;
    clang_visitChildren(
        c,
        [](CXCursor c, CXCursor _, CXClientData cd) {
          switch (clang_getCursorKind(c)) {
          case CXCursor_ParmDecl: {

            vector<string> *ret = (vector<string> *)cd;

            string parmtypestring =
                typeToStringOCaml(clang_getCursorType(c), c);
            ret->push_back(parmtypestring);
            break;
          }
          default:
            break;
          }
          return CXChildVisit_Continue;
        },
        (CXClientData)&parmtypes);
    string res = "Ctypes.(";
    for (auto pt : parmtypes) {
      res += (pt + " @-> ");
    }
    res += "returning (" + rettype + ")) static_funptr";
    return res;
  }
  case CXType_Unexposed: {
    // Function pointers inside a struct result in this case. Getting
    // the canonical type returns a function proto type which is
    // processed in the recursive call.
    CXType t = clang_getCanonicalType(ty);
    return typeToStringOCaml(t, c);
  }
  case CXType_Pointer: {
    CXType ptee = clang_getPointeeType(ty);
    switch (ptee.kind) {
    case CXType_FunctionProto:
      return "(" + typeToStringOCaml(ptee, c) + ")";
    case CXType_Char_U:
    case CXType_Char_S:
      return "string";
    default:
      return "(" + typeToStringOCaml(ptee, c) + ") ptr";
    }
  }
  default: {
    string s = getFromCXString(clang_getTypeSpelling(ty));
    replaceAll(s, std::string("const struct"), std::string(""));
    replaceAll(s, std::string("struct"), std::string(""));
    return s;
  }
  }
}

string typeToString(CXType ty, CXCursor c) {
  switch (ty.kind) {
  case CXType_Void:
    return "void";
  case CXType_UInt:
    return "size_t";
  case CXType_Int:
    return "int";
  case CXType_Bool:
    return "bool";
  case CXType_ULong:
    return "ulong";
  case CXType_UChar:
    return "uchar";
  case CXType_Char_S:
  case CXType_Char_U:
    return "char";
  case CXType_ConstantArray: {
    string s = getFromCXString(clang_getTypeSpelling(ty));
    return "array " + to_string(clang_getArraySize(ty)) + " (" +
           typeToString(clang_getArrayElementType(ty), c) + ")";
  }
  case CXType_IncompleteArray: {
    string s = getFromCXString(clang_getTypeSpelling(ty));
    return "ptr (" + typeToString(clang_getArrayElementType(ty), c) + ")";
  }
  case CXType_Typedef: {
    string s = getFromCXString(clang_getTypeSpelling(ty));
    // Ugly hack: check if the typedef contains an int and if so, use
    // the OCaml equivalent.
    size_t found = s.find("int");
    if (found != std::string::npos) {
      return s;
    }
    // Otherwise return the typedef'd type.
    return getFromCXString(clang_getTypeSpelling(ty));
  }
  case CXType_FunctionProto: {
    // For function proto we need to do some extra work by traversing
    // the cursor to get the parameter types.

    string rettype = typeToString(clang_getResultType(ty), c);
    vector<string> parmtypes;
    clang_visitChildren(
        c,
        [](CXCursor c, CXCursor p, CXClientData cd) {
          switch (clang_getCursorKind(c)) {
          case CXCursor_ParmDecl: {

            vector<string> *ret = (vector<string> *)cd;

            string parmtypestring = typeToString(clang_getCursorType(c), c);
            ret->push_back(parmtypestring);
            break;
          }
          default:
            break;
          }
          return CXChildVisit_Continue;
        },
        (CXClientData)&parmtypes);
    string res = "static_funptr Ctypes.(";
    for (auto pt : parmtypes) {
      res += (pt + " @-> ");
    }
    res += "returning (" + rettype + "))";
    return res;
  }
  case CXType_Unexposed: {
    // Function pointers inside a struct result in this case. Getting
    // the canonical type returns a function proto type which is
    // processed in the recursive call.
    CXType t = clang_getCanonicalType(ty);
    return typeToString(t, c);
  }
  case CXType_Pointer: {
    CXType ptee = clang_getPointeeType(ty);
    switch (ptee.kind) {
    case CXType_FunctionProto: {
      return "(" + typeToString(ptee, c) + ")";
    }
    case CXType_Char_U:
    case CXType_Char_S: {
      return "string";
    }
    default:
      return "ptr (" + typeToString(ptee, c) + ")";
    }
  }
  default: {
    string s = getFromCXString(clang_getTypeSpelling(ty));
    replaceAll(s, std::string("const struct"), std::string(""));
    replaceAll(s, std::string("struct"), std::string(""));
    return s;
  }
  }
}

CXChildVisitResult gatherStructDecls(CXCursor c, CXCursor parent,
                                     CXClientData client_data) {
  if (clang_Location_isInSystemHeader(clang_getCursorLocation(c)) != 0) {
    return CXChildVisit_Continue;
  }
  switch (clang_getCursorKind(c)) {
  case CXCursor_TypedefDecl: {
    if (clang_getTypedefDeclUnderlyingType(c).kind != CXType_Elaborated) {
      break;
    }
    // fall through (typedef'd type can be a struct)
  }
  case CXCursor_StructDecl: {
    vector<StructDecl> *v = (vector<StructDecl> *)client_data;
    StructDecl s;
    s.name = getFromCXString(clang_getCursorSpelling(c));
    auto already_visited = [&](StructDecl d) {
      return d.name.compare(s.name) == 0;
    };
    auto found = std::find_if(std::begin(*v), std::end(*v), already_visited);
    if (found != std::end(*v)) {
      break;
    }
    (clang_getCursorKind(c) == CXCursor_TypedefDecl) ? s.is_typedef = true
                                                     : s.is_typedef = false;
    tuple<string, vector<tuple<string, string>> *> data = {s.name,
                                                           &s.fieldnames};
    clang_visitChildren(
        c,
        [](CXCursor c, CXCursor p, void *d) {
          switch (clang_getCursorKind(c)) {
          case CXCursor_FieldDecl: {
            CXType f_type = clang_getCursorType(c);
            tuple<string, vector<tuple<string, string>> *> *data =
                (tuple<string, vector<tuple<string, string>> *> *)d;
            vector<tuple<string, string>> *fieldnames = std::get<1>(*data);
            string sname = std::get<0>(*data);
            string fieldname = (getFromCXString(clang_getCursorSpelling(c)));
            string type = lowerCase(typeToString(f_type, c));
            replaceAll(type, std::string("ctypes.("), std::string("T.("));
            fieldnames->push_back(make_tuple(fieldname, type));
            break;
          }
          default:
            break;
          }
          return CXChildVisit_Recurse;
        },
        (CXClientData)&data);
    v->push_back(s);
    break;
  }
  default:
    break;
  }
  return CXChildVisit_Recurse;
}

CXChildVisitResult gatherFuncDecls(CXCursor c, CXCursor parent,
                                   CXClientData client_data) {
  if (clang_Location_isInSystemHeader(clang_getCursorLocation(c)) != 0) {
    return CXChildVisit_Continue;
  }
  switch (clang_getCursorKind(c)) {
  case CXCursor_FunctionDecl: {
    vector<FuncDecl> *v = (vector<FuncDecl> *)client_data;
    FuncDecl f;
    CXType fty = clang_getCursorType(c);
    CXType retty = clang_getResultType(fty);
    f.name = getFromCXString(clang_getCursorSpelling(c));
    f.resulttype = lowerCase(typeToString(retty, c));
    replaceAll(f.resulttype, std::string("ctypes.("), std::string("Ctypes.("));

    tuple<string, vector<string> *> data = {f.name, &f.paramtypes};
    clang_visitChildren(
        c,
        [](CXCursor c, CXCursor p, CXClientData d) {
          tuple<string, vector<string> *> *f =
              (tuple<string, vector<string> *> *)d;
          vector<string> *paramtypes = get<1>(*f);
          string name = get<0>(*f);
          switch (clang_getCursorKind(c)) {
          case CXCursor_ParmDecl: {
            CXType c_type = clang_getCursorType(c);
            CXCursorKind p_type = clang_getCursorKind(p);
            string s = lowerCase(typeToString(c_type, c));
            replaceAll(s, std::string("ctypes.("), std::string("Ctypes.("));
            // To avoid parsing the function pointer
            if (p_type != CXCursor_ParmDecl) {
              paramtypes->push_back(s);
            }
            break;
          }
          default:
            break;
          }
          return CXChildVisit_Recurse;
        },
        (CXClientData)&data);
    if (f.paramtypes.size() == 0) {
      f.paramtypes.push_back("void");
    }
    v->push_back(f);
    break;
  }
  default:
    break;
  }
  return CXChildVisit_Recurse;
}

CXChildVisitResult gatherTypedefDecls(CXCursor c, CXCursor parent,
                                      CXClientData client_data) {
  if (clang_Location_isInSystemHeader(clang_getCursorLocation(c)) != 0) {
    return CXChildVisit_Continue;
  }
  switch (clang_getCursorKind(c)) {
  case CXCursor_TypedefDecl: {
    vector<TypedefDecl> *v = (vector<TypedefDecl> *)client_data;
    TypedefDecl t;
    t.name = getFromCXString(clang_getCursorSpelling(c));
    auto is_type_alias = true;
    auto typ = clang_getTypedefDeclUnderlyingType(c);
    string s = (typeToString(typ, c));
    t.type.append(s);
    t.ocaml_type = typeToStringOCaml(typ, c);
    clang_visitChildren(
        c,
        [](CXCursor c, CXCursor p, void *d) {
          bool *is_type_alias = (bool *)d;
          switch (clang_getCursorKind(c)) {
          case CXCursor_StructDecl:
          case CXCursor_FunctionDecl:
          case CXCursor_EnumConstantDecl:
            *is_type_alias = false;
            break;
          default:
            break;
          }
          return CXChildVisit_Recurse;
        },
        (CXClientData)&is_type_alias);

    if (is_type_alias) {
      v->push_back(t);
    }
    break;
  }
  default:
    break;
  }
  return CXChildVisit_Recurse;
};

CXChildVisitResult gatherEnumDecl(CXCursor c, CXCursor parent,
                                  CXClientData client_data) {
  if (clang_Location_isInSystemHeader(clang_getCursorLocation(c)) != 0) {
    return CXChildVisit_Continue;
  }
  switch (clang_getCursorKind(c)) {
  case CXCursor_EnumDecl: {
    vector<EnumDecl> *v = (vector<EnumDecl> *)client_data;
    EnumDecl e;
    e.name = getFromCXString(clang_getCursorSpelling(c));
    clang_visitChildren(
        c,
        [](CXCursor c, CXCursor p, CXClientData cd) {
          switch (clang_getCursorKind(c)) {
          case CXCursor_EnumConstantDecl: {
            vector<string> *e = (vector<string> *)cd;
            string s = getFromCXString(clang_getCursorSpelling(c));
            e->push_back(s);
            break;
          }
          default:
            break;
          }
          return CXChildVisit_Recurse;
        },
        (CXClientData)&e.enumconstants);
    v->push_back(e);
    break;
  }
  default:
    break;
  }
  return CXChildVisit_Recurse;
};

vector<EnumDecl> generateEnumBindings(CXCursor cur) {
  vector<EnumDecl> enums;
  clang_visitChildren(cur, gatherEnumDecl, &enums);
  return enums;
}

vector<StructDecl> generateStructBindings(CXCursor cur) {
  vector<StructDecl> structdecls;
  clang_visitChildren(cur, gatherStructDecls, &structdecls);
  return structdecls;
}

vector<FuncDecl> generateFuncBindings(CXCursor cur) {
  vector<FuncDecl> funcdecls;
  clang_visitChildren(cur, gatherFuncDecls, &funcdecls);
  return funcdecls;
}

vector<TypedefDecl> generateTypedefBindings(CXCursor cur) {
  vector<TypedefDecl> typedefdecls;
  clang_visitChildren(cur, gatherTypedefDecls, &typedefdecls);
  return typedefdecls;
}

template <typename T>
void removeDecls(vector<string> blacklist, vector<T> *decls) {
  sort(begin(blacklist), end(blacklist));
  decls->erase(remove_if(begin(*decls), end(*decls),
                         [&](auto x) {
                           return binary_search(begin(blacklist),
                                                end(blacklist), x);
                         }),
               end(*decls));
}

void generateBindings(CXCursor cur, std::ofstream &types_out,
                      std::ofstream &funcs_out) {
  vector<StructDecl> strucs = generateStructBindings(cur);
  vector<FuncDecl> funcs = generateFuncBindings(cur);
  vector<EnumDecl> enums = generateEnumBindings(cur);
  vector<TypedefDecl> typedefs = generateTypedefBindings(cur);

  removeDecls(blacklist, &funcs);
  removeDecls(blacklist, &strucs);

  types_out << "open Ctypes\n";

  types_out << "module Types (T : Ctypes.TYPE) ="
            << "\n";
  types_out << "struct"
            << "\n";
  types_out << "\t"
            << "open T"
            << "\n";

  // print out type definitions first
  for (auto t : typedefs) {

    string typedef_lowername = lowerCase(t.name);

    types_out << "\ttype " << typedef_lowername << " = " << t.ocaml_type
              << "\n";
    replaceAll(t.type, std::string("Ctypes.("), std::string("T.("));

    types_out << "\tlet " << typedef_lowername << " : " << typedef_lowername
              << " typ = " << t.type << "\n";
    types_out << "\n";
  }

  for (auto e : enums) {
    if (e.name.empty()) {
      continue;
    }

    string enum_lowername = lowerCase(e.name);

    types_out << "\ttype " << enum_lowername << " ="
              << "\n";
    for (auto ec : e.enumconstants) {
      types_out << "\t\t | " << upperCase(ec) << "\n";
    }
    types_out << "\n";
  }

  for (auto s : strucs) {
    if (s.name.empty()) {
      continue;
    }

    types_out << "\ttype " << lowerCase(s.name) << "\n";
  }
  types_out << "\n";

  for (auto e : enums) {
    if (e.name.empty()) {
      continue;
    }

    string enum_lowername = lowerCase(e.name);
    for (auto s : e.enumconstants) {
      types_out << "\t"
                << "let " << lowerCase(s) << " = constant \"" << s
                << "\" int64_t"
                << "\n";
    }
    types_out << endl;

    types_out << "\t"
              << "let " << enum_lowername << " = enum \"" << e.name << "\" ["
              << "\n";

    for (auto s : e.enumconstants) {
      types_out << "\t\t(" << upperCase(s) << ", " << lowerCase(s) << ");"
                << "\n";
    }
    types_out << "\t"
              << "]" << endl
              << "\n";
  }

  for (auto s : strucs) {
    if (s.name.empty()) {
      continue;
    }
    string struct_lowername = lowerCase(s.name);

    if (s.is_typedef) {
      types_out << "\tlet " << struct_lowername << " : " << struct_lowername
                << " structure typ = typedef (structure \"" << s.name << "\")"
                << " \"" << s.name << "\""
                << "\n";
    } else {
      types_out << "\tlet " << struct_lowername << " : " << struct_lowername
                << " structure typ = structure \"" << s.name << "\""
                << "\n";
    }
    for (auto f : s.fieldnames) {
      replaceAll(get<1>(f), std::string("Ctypes.("), std::string("T.("));

      types_out << "\tlet " << struct_lowername + "_" + get<0>(f) << " = field "
                << struct_lowername << " \"" << get<0>(f) << "\" (" << get<1>(f)
                << ")"
                << "\n";
    }
    types_out << "\tlet () = seal " << struct_lowername << "\n\n";
  }

  types_out << "end";

  funcs_out << "open Ctypes\n"
            << "open Types_generated\n";

  int count = 0;
  int n = 0;

  for (auto f : funcs) {

    if (count % 100 == 0) {
      if (count == 0) {
        funcs_out << "module F0 "
                     "(F : Ctypes.FOREIGN) ="
                  << "\n";
        funcs_out << "struct"
                  << "\n";
        funcs_out << "\topen F"
                  << "\n";

      } else {

        funcs_out << "end\n";
        n++;

        funcs_out << "module F" << n << " (F : Ctypes.FOREIGN) ="
                  << "\n";
        funcs_out << "struct"
                  << "\n";
        funcs_out << "\topen F"
                  << "\n";
      }
    }

    funcs_out << "\tlet " << lowerCase(f.name) << " = foreign \"" << f.name
              << "\" (";
    for (auto e : f.paramtypes) {
      funcs_out << e << " @-> ";
    }

    funcs_out << "returning (" << f.resulttype << "))"
              << "\n";
    count++;
  }
  funcs_out << "end\n";

  funcs_out << "module Functions (F : Ctypes.FOREIGN) = struct\n";
  for (int i = 0; i < n; i++) {
    funcs_out << "\tinclude F" << i << " (F)\n";
  }
  funcs_out << "end";
}

void configureParser(cli::Parser &parser) {
  parser.set_required<std::string>("f", "file", "", "header file to process");
  parser.set_optional<bool>("c++", "c++", false, "Parse file as a C++ file");
}

// https://stackoverflow.com/questions/62005698/libclang-clang-getargtype-returns-wrong-type/66702426#66702426
void printDiagnostics(CXTranslationUnit translationUnit) {
  int nbDiag = clang_getNumDiagnostics(translationUnit);
  printf("There are %i diagnostics:\n", nbDiag);

  bool foundError = false;
  for (unsigned int currentDiag = 0; currentDiag < nbDiag; ++currentDiag) {
    CXDiagnostic diagnotic = clang_getDiagnostic(translationUnit, currentDiag);
    CXString errorString = clang_formatDiagnostic(
        diagnotic, clang_defaultDiagnosticDisplayOptions());
    std::string tmp{clang_getCString(errorString)};
    clang_disposeString(errorString);
    if (tmp.find("error:") != std::string::npos) {
      foundError = true;
    }
    std::cerr << tmp << std::endl;
  }
}

int main(int argc, char *argv[]) {
  cli::Parser parser(argc, argv);
  configureParser(parser);
  parser.run_and_exit_if_error();

  string filename = parser.get<string>("f");

  CXIndex index = clang_createIndex(0, 0);
  auto unit = parseFile(index, filename, parser.get<bool>("c++"));

  CXCursor cur = clang_getTranslationUnitCursor(unit);
  std::ofstream types_out;
  std::ofstream funcs_out;
  types_out.open("type_description.ml", ios::out | ios::trunc);
  funcs_out.open("function_description.ml");
  generateBindings(cur, types_out, funcs_out);
  types_out.close();
  funcs_out.close();

  printDiagnostics(unit);

  clang_disposeTranslationUnit(unit);
  clang_disposeIndex(index);

  return 0;
}
