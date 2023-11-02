#pragma once
#include "THStack.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveText.h"

std::vector<int> colors = {kRed, kBlue, kGreen, kBlack, kYellow, kMagenta, kCyan};

template<typename T> class StackWrapper{
private:
	THStack* hs;
	std::vector<TH1F*> hists;
	TLegend* leg;
	TString name;
	std::string title;

	uint nbins;
	float xlow, xhigh;
	bool fill;
public:
	StackWrapper()
	: hs(nullptr), name(""), title(""), nbins(0), xlow(0), xhigh(0), fill(false){};
	StackWrapper(std::string _name, std::string _title, uint _nbins, float _xlow, float _xhigh, bool _fill = true,
		bool legend_left = false, Style_t linestyle = kSolid, Style_t fillstyle=kFEmpty)
		: name(_name), title(_title), nbins(_nbins), xlow(_xlow), xhigh(_xhigh), fill(_fill){

		std::string hstitle = title.substr(0, title.find(';'));

		hs = new THStack(hstitle.c_str(), title.c_str());

		if(legend_left)
			leg = new TLegend(0.1,0.7,0.22,0.9);
		else
			leg = new TLegend(0.58,0.7,0.9,0.9);


		for(uint p = 0; p < T::UNDEF; p++){
			TH1F* h = new TH1F("", "", nbins, xlow, xhigh);
			if(fill){
				h->SetFillColor(colors[p]);
				h->SetLineColorAlpha(colors[p], 0.5);
			}
			else
				h->SetLineColor(colors[p]);
			h->SetLineStyle(linestyle);
			h->SetFillStyle(fillstyle);
			hists.push_back(h);
			hs->Add(h);

			leg->AddEntry(h, T::names[p], "f");
		}
	};
	virtual ~StackWrapper(){
		// for(TH1F *h : hists)
		// 	delete h;
		// delete hs;
		// delete leg;
	};
	THStack* Stack(){return hs;};
	void Fill(int p, Double_t value){
		if(p < T::UNDEF){
			hists[p]->Fill(value);
		}
	};

	void Fill(int p, Double_t value, Double_t weight){
		if(p < T::UNDEF){
			hists[p]->Fill(value, weight);
		}
	};

	TCanvas* Draw(TString draw_option = "", TCanvas *c = nullptr, bool addTotal = false){
		if(c == nullptr){
			c = new TCanvas();
		}
		c->cd();
		hs->Draw(draw_option);

		if(addTotal){
			TH1* sum = static_cast<TH1*>(hs->GetStack()->Last());
			sum->SetLineStyle(kDashed);
			sum->SetLineColor(kBlack);
			sum->Draw("SAME");
			leg->AddEntry(sum, "Sum");
		}

		leg->Draw();
		return c;
	};

	void Write(){
		hs->Write();
	};

	void normalize(){
		for(TH1F* h : hists){
			if(h->GetEntries() != 0)
				h->Scale(1./h->GetEntries());
		}
	};

	void SetTitle(TString new_title){
		title = new_title;
		hs->SetTitle(new_title);
	}

	TLegend* Legend(){
		return leg;
	}

	std::vector<TH1F*> GetHists(){
		return hists;
	}

	TH1F* Get(int p){
		return hists[p];
	}

	TString GetName() const{
		return name;
	}

	TPaveText* GetStats(){
		TPaveText *p = new TPaveText(0.1,0.7,0.22,0.9, "NDC");
		uint i = 0;
		for(TH1F* h : hists){
			float mean = h->GetMean();
			float std = h->GetStdDev();

			TText *text = p->AddText(Form("#mu = %.2f ; #sigma = %.2f", mean, std));
			text->SetTextColor(colors[i]);
			i++;
		}
		return p;
	}
};