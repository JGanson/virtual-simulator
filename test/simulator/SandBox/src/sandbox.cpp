#include "Header.h"

class Sandbox : public temp::Application {
public:
	Sandbox() {

	}

	~Sandbox() {
	
	}
};

temp::Application*	temp::Creat_application() {
	return new Sandbox();
}