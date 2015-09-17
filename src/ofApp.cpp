#include "ofApp.h"

//--------------------------------------------------------------
void ofApp::setup(){

    D_alpha = 1.0;
    alpha = 0.75;

    dX = 1.0;
    dT = powf(dX * dX / (2.0 * D_alpha), 1.0 / alpha); // Calc here!
   
    base_dT = dT;

    // Now omega is calculated to give the same D_alpha 
    omega = 2.0 * dT * D_alpha / (dX * dX);
    r = 1.0;
    
    max_n = MAX_N;
    survival_sib = new double[max_n];
    survival_exp = new double[max_n];
    survival_sib[0] = 1.0 - alpha;
    survival_exp[0] = 1.0 - omega;
    for (int i=1; i<max_n; i++) {
        survival_sib[i] = survival_sib[i-1] * (1. - alpha/((double)(i+1)));;
        survival_exp[i] = survival_exp[i-1] * (1.0 * omega);
    }
    sib_stats = new double[max_n];
    exp_stats = new double[max_n];
    sib_dist = new double[max_n];
    sib_dist[0] = alpha;
    for (int i=1; i<max_n; i++) 
        sib_dist[i] = (double(i)/double(i+1)) * (1.0 - alpha / double(i)) * sib_dist[i-1];
    sjc = 0;

    n = 0;
    t = 0.0;
    sib_pos = ofPoint::zero(); //= {0, 0, 0};
    exp_pos = ofPoint::zero(); //{0, 0, 0};

    sib_pos_hist.push_back(sib_pos);
    exp_pos_hist.push_back(exp_pos);
    
    seed = -1;
    double sib_rand = ran2(&seed);
    double exp_rand = ran2(&seed);
    sib_arrival = inv_cum_dist(sib_rand, survival_sib, max_n);
    exp_arrival = inv_cum_dist(exp_rand, survival_exp, max_n);
    sib_arrival_t = sib_arrival * dT;
    exp_arrival_t = exp_arrival * dT;
    
    sib_color = ofColor(0,255,0,200);//::green;
    exp_color = ofColor(255,0,0,200);//::red;
    
    /*
    sib_graph.setMode(OF_PRIMITIVE_LINE_STRIP);
    sib_graph.addColor(sib_color);
    sib_graph.addVertex(ofVec2f(0.5 * ofGetWidth() + sib_pos[0] * dX,
                                0.5 * ofGetHeight() + sib_pos[1] * dX));

    exp_graph.setMode(OF_PRIMITIVE_LINE_STRIP);
    exp_graph.addColor(exp_color);
    exp_graph.addVertex(ofVec2f(0.5 * ofGetWidth() + exp_pos[0] * dX,
                                0.5 * ofGetHeight() + exp_pos[1] * dX));
    */

    ofSetWindowTitle("Comparison of anomalous and regular diffusion");
    ofDisableAntiAliasing();
    ofSetFrameRate(FRAME_RATE);

	gui.setup(); // most of the time you don't need a name
	gui.add(rate.setup("time stretch", FRAME_RATE, 1, 200 ));
    gui.setPosition(10, ofGetHeight() - gui.getHeight() - 10);
}

//--------------------------------------------------------------
void ofApp::update(){
    
    ofSetFrameRate(rate);

    n++;
    t += 100.0 / float(rate);

    while (n >= exp_arrival) {
        exp_arrival += exp_wait;
        exp_arrival_t += exp_wait * dT;

        //exp_graph.addColor(exp_color);
        rand_direction(ran2(&seed), ran2(&seed), exp_pos);
        exp_pos_hist.push_back(exp_pos);

        double exp_rand = ran2(&seed);
        //exp_arrival += inv_cum_dist(exp_rand, survival_exp, max_n);
        exp_wait = inv_cum_dist(exp_rand, survival_exp, max_n);
        
        exp_cd_dt = dT;

        //exp_graph.addVertex(ofVec2f(0.5 * ofGetWidth() + exp_pos[0] * dX,
        //                            0.5 * ofGetHeight() + exp_pos[1] * dX));
    }
    while (n >= sib_arrival) {
        sib_arrival += sib_wait;
        sib_arrival_t += sib_wait * dT;

        rand_direction(ran2(&seed), ran2(&seed), sib_pos);
        sib_pos_hist.push_back(sib_pos);
        
        double sib_rand = ran2(&seed);
        sib_wait = inv_cum_dist(sib_rand, survival_sib, max_n);
        
        sib_cd_dt = dT;

        sib_stats[sib_wait-1] += 1.0;
        sjc++;

        //sib_graph.addColor(sib_color);
        //sib_graph.addVertex(ofVec2f(0.5 * ofGetWidth() + sib_pos[0] * dX,
        //                            0.5 * ofGetHeight() + sib_pos[1] * dX));
    }


    
}

//--------------------------------------------------------------
void ofApp::draw(){
	ofBackgroundGradient(ofColor(50), ofColor(0));
    

    ofPoint midpoint(0.5 * ofGetWidth(), 0.5 * ofGetHeight(), 0);
    
    ofMesh sib_graph;
    sib_graph.setMode(OF_PRIMITIVE_LINE_STRIP);
    for (int i=0; i<sib_pos_hist.size(); i++) {
        sib_graph.addColor(sib_color);
        sib_graph.addVertex(midpoint + sib_pos_hist[i] * dX);
    }                

    ofMesh exp_graph;
    exp_graph.setMode(OF_PRIMITIVE_LINE_STRIP);
    for (int i=0; i<exp_pos_hist.size(); i++) {
        exp_graph.addColor(exp_color);
        exp_graph.addVertex(midpoint + exp_pos_hist[i] * dX);
    }   
    
    exp_graph.draw();
    // Little marker circle
    ofSetColor(exp_color);
    ofNoFill();
    ofCircle(0.5 * ofGetWidth() + exp_pos[0] * dX,
             0.5 * ofGetHeight() + exp_pos[1] * dX, 5);
    ofSetColor(ofColor(255,255,255,255));
    ofFill();
    ofCircle(0.5 * ofGetWidth() + exp_pos[0] * dX,
             0.5 * ofGetHeight() + exp_pos[1] * dX, 4);
    
    ofSetColor(ofColor(255));
    ofDrawBitmapString(ofToString(exp_arrival + exp_wait - n - 1), midpoint[0] + exp_pos[0] * dX + 6, midpoint[1] + exp_pos[1] * dX - 3);

    sib_graph.draw();
    // Little marker circle
    ofSetColor(sib_color);
    ofNoFill();
    ofCircle(0.5 * ofGetWidth() + sib_pos[0] * dX,
             0.5 * ofGetHeight() + sib_pos[1] * dX, 5);
    ofSetColor(ofColor(255,255,255,255));
    ofFill();
    ofCircle(0.5 * ofGetWidth() + sib_pos[0] * dX,
             0.5 * ofGetHeight() + sib_pos[1] * dX, 4);

    ofSetColor(ofColor(255));
    ofDrawBitmapString(ofToString(sib_arrival + sib_wait - n - 1), midpoint[0] + sib_pos[0] * dX + 6, midpoint[1] + sib_pos[1] * dX - 3);
 
    ofSetColor(255);
    ofCircle(0.5 * ofGetWidth(), 0.5 * ofGetHeight(), 5);


    ofSetColor(ofColor(255));
    ofDrawBitmapString("Step " + ofToString(n) + ", alpha=" + ofToString(alpha) 
                       + " dT=" + ofToString(dT) + ", dX=" + ofToString(dX) + " omega="
                       + ofToString(omega), 
                       TW_MARGIN, TH_MARGIN + 1*12);
    ofSetColor(sib_color);
    ofDrawBitmapString("Anomalous Diffusion Position [" + ofToString(sib_pos[0]) + ", "
                       + ofToString(sib_pos[1]) + ", "
                       + ofToString(sib_pos[2]) + "]",
                       TW_MARGIN, TH_MARGIN + 2*12);
    ofSetColor(exp_color);
    ofDrawBitmapString("Brownian Motion Position [" + ofToString(exp_pos[0]) + ", "
                       + ofToString(exp_pos[1]) + ", "
                       + ofToString(exp_pos[2]) + "]",
                       TW_MARGIN, TH_MARGIN + 3*12);

    /*ofSetColor(ofColor(255));
    ofDrawBitmapString(ofToString(sib_dist[0]) +", "+ ofToString(sib_dist[1]) +", "+ ofToString(sib_dist[2]) +", "+  
                       ofToString(sib_dist[3]) +", "+ ofToString(sib_dist[4]) +", "+ ofToString(sib_dist[5]), 
                       TW_MARGIN, TH_MARGIN + 4*12);
    ofDrawBitmapString(ofToString(sib_stats[0]/double(sjc)) +", "+ ofToString(sib_stats[1]/double(sjc)) +", "+ ofToString(sib_stats[2]/double(sjc)) +", "+  
                       ofToString(sib_stats[3]/double(sjc)) +", "+ ofToString(sib_stats[4]/double(sjc)) +", "+ ofToString(sib_stats[5]/double(sjc)), 
                       TW_MARGIN, TH_MARGIN + 5*12);
*/
    gui.draw();
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){

}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){
    dX *= (1.0 + 2.0 * double(x - mouse_x_start) / double(ofGetHeight()));
    dX = max(dX, 0.4);
    dT = powf(dX * dX / (2.0 * D_alpha), 1.0 / alpha); // Calc here!
}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){
    mouse_x_start = x;
}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

}
