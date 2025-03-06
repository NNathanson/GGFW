#include "../BootstrapProfile.h"
#include <iostream>

  void test(){ 
    gSystem->Load("../libGGFW");
    BootstrapProfile* bs = new BootstrapProfile;
    std::cout<<"testing"<<std::endl;
    return;
  }
