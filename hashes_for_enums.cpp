#include <unordered_map>
#include <tuple>
#include <iostream>

enum VARIABLE_NATURE{
NEW_SOLIDS,
OLD_SOLIDS,
BLANK,
ERROR
};

enum VARIABLE_ROLE{
OWNED,
OVERLAP
};
/*

namespace std {
	template<class E>class hash {
		using bogotype = typename std::enable_if<std::is_enum<E>::value, E>::type;
		public:
		size_t operator()(const E&e) const {
			return std::hash<typename std::underlying_type<E>::type>()(e);
		}
	};

}
*/
namespace std {
	//Full specialization of hash for my enums
	template<> struct hash<VARIABLE_NATURE> {
	 using bogotype = typename std::enable_if<std::is_enum<VARIABLE_NATURE>::value, VARIABLE_NATURE>::type;
			  public:
		    size_t operator()(const VARIABLE_NATURE&e) const {
			return std::hash<typename std::underlying_type<VARIABLE_NATURE>::type>()(e);
			  }
	 };

	template<> struct hash<VARIABLE_ROLE> {
	 using bogotype = typename std::enable_if<std::is_enum<VARIABLE_ROLE>::value, VARIABLE_ROLE>::type;
			  public:
		    size_t operator()(const VARIABLE_ROLE&e) const {
			return std::hash<typename std::underlying_type<VARIABLE_ROLE>::type>()(e);
			  }
	 };

}



int main(int argc, char **argv){
	std::unordered_map<std::string, std::pair<VARIABLE_ROLE, VARIABLE_NATURE> > varNameToVarArchetypeDict; // keys are variable names
	varNameToVarArchetypeDict["apple_vector_extreme"] = std::pair<VARIABLE_ROLE, VARIABLE_NATURE>(OVERLAP, NEW_SOLIDS);

	std::unordered_map<VARIABLE_ROLE, std::string> roleEnumToRoleNameDict = {{OWNED,"OWNED"}, {OVERLAP,"OVERLAP"}};
	std::unordered_map<VARIABLE_NATURE, std::string> natureEnumToNatureNameDict = {{NEW_SOLIDS,"NEW_SOLIDS"}, {OLD_SOLIDS, "OLD_SOLIDS"}, {BLANK, "BLANK"}, {ERROR,"ERROR"}};

	std::cout << "The role of apple_vector_extreme is: " << roleEnumToRoleNameDict[varNameToVarArchetypeDict["apple_vector_extreme"].first] <<
		  " and the nature of apple_vector_extreme is: " << natureEnumToNatureNameDict[varNameToVarArchetypeDict["apple_vector_extreme"].second] << std::endl; 

	return 0;	
}
