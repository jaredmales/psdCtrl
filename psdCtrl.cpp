
#include <mx/app/application.hpp>
using namespace mx::app;

#include <mx/ioutils/readColumns.hpp>
using namespace mx::ioutils;

#include <mx/sigproc/averagePeriodogram.hpp>
#include <mx/sigproc/signalWindows.hpp>
using namespace mx::sigproc;

#include <mx/math/vectorUtils.hpp>
using namespace mx::math;

template<typename _realT>
class psdCtrl : public application
{

public:

    typedef _realT realT;

public:

    std::string m_infile;

    std::string m_sep {"space"};

    std::string m_outfile {"psd_out.dat"};

    std::string m_windowName {"hann"};

    std::vector<realT> m_time;

    std::vector<realT> m_dist;

    virtual void setupConfig();

    virtual void loadConfig();

    virtual int execute();
    
};

template<typename realT>
void psdCtrl<realT>::setupConfig()
{
    config.add("input.file", "i", "input.file", argType::Required, "input", "file", false, "string", "input file.");
    config.add("input.sep", "", "input.sep", argType::Required, "input", "sep", false, "string", "column separator, either 'space' (default) or 'comma'.");
    config.add("output.file", "o", "output.file", argType::Required, "output", "file", false, "string", "output file.");
    config.add("window.name", "w", "window.name", argType::Required, "window", "name", false, "string", "name of the window.  Possible values: hann (default), rect");

    m_nonOptionHelp = "[input-file] [output-file]";
}

template<typename realT>
void psdCtrl<realT>::loadConfig()
{
    config(m_infile, "input.file");

    if(config.nonOptions.size() > 0)
    {
        m_infile = config.nonOptions[0];
    }

    config(m_sep, "input.sep");

    config(m_outfile, "output.file");
    
    if(config.nonOptions.size() > 1)
    {
        m_outfile = config.nonOptions[1];
    }

    config(m_windowName, "window.name");
}

template<typename realT>
int psdCtrl<realT>::execute()
{
    if(m_infile == "")
    {
        std::cerr << "no input file provided\n";
        help();
        return -1;
    }

    m_time.clear();
    m_dist.clear();
    if(m_sep == "space")
    {
        if(readColumns(m_infile, m_time, m_dist) < 0)
        {
            std::cerr << "Error reading file\n";
            help();
            return -1;
        }
    }
    else if(m_sep == "comma")
    {
        if(readColumns<','>(m_infile, m_time, m_dist) < 0)
        {
            std::cerr << "Error reading file\n";
            help();
            return -1;
        }

    }
    else
    {
        std::cerr << "invalid separator specified.  must be either 'space' or 'comma'\n";
        help();
        return -1;
    }
    
    std::cout << "rms: " << sqrt(vectorVariance(m_dist)) << "\n";

    averagePeriodogram<realT> avper(m_time.size(), 0, m_time[1]-m_time[0]);

    if(m_windowName == "hann")
    {
        avper.win(window::hann);
    }

    std::vector<realT> psd;

    avper(psd, m_dist);

    std::vector<realT> freq(psd.size());
    for(size_t n = 0; n < psd.size(); ++n)
    {
        freq[n] = avper[n];
    }

    realT m_maxFreq = 2000;
    realT m_freqAlpha = -3;

    if(freq.back() < m_maxFreq)
    {
        realT df = freq[1] - freq[0];

        while(freq.back() < m_maxFreq)
        {
            freq.push_back(freq.back() + df);
            psd.push_back( psd.back() * pow( freq.back()/(freq.back()-df), m_freqAlpha));
        }
    }



    std::ofstream fout;
    fout.open(m_outfile);
    for(size_t n = 0; n < psd.size(); ++n)
    {
        fout << avper[n] << " " << psd[n] << "\n";
    }
    fout.close();

    return 0;

}

int main( int argc,
          char ** argv
        )
{
    psdCtrl<float> psdc;

    return psdc.main(argc, argv);
}
