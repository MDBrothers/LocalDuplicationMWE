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

		namespace std {
	 		template<class E>class hash {
			 using  = typename std::enable_if<std::is_enum<E>::value, E>::type;
					  public:
				    size_t operator()(const E&e) const {
				    	return std::hash<typename std::underlying_type<E>::type>()(e);
					  }
			 };
		}

		std::unordered_map<VARIABLE_ROLE, VARIABLE_NATURE > varNameToVarArchetypeDict; // keys are variable names

