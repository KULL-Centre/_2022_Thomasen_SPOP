import hoomd
import hoomd.md
import time
import os
import sys
import itertools
import pandas as pd
import numpy as np
import mdtraj as md
from hoomd import azplugins
from argparse import ArgumentParser
from scipy.optimize import curve_fit
from itertools import combinations
sys.path.append('/groups/sbinlab/thomasen/software/BLOCKING')
from main import BlockAnalysis

parser = ArgumentParser()
parser.add_argument('--seq_name',nargs='?',required=True)
parser.add_argument('--path',nargs='?',required=True)
args = parser.parse_args()

print(hoomd.__file__)

def initProteins():
    proteins = pd.DataFrame(index=['SETD2_IDR1', 'SETD2_IDR2', 'SCAF1_IDR1', 'SRC3_IDR1', 'Gli3_IDR1', 'Gli3_IDR2', 'Gli2_IDR2'], columns=['eps_factor','pH','ionic','fasta'])
    fasta_SETD2_IDR1 = """MKQLQPQPPPKMGDFYDPEHPTPEEEENEAKIENVQKTGFIKGPMFKGVASSRFLPKGTKTKVNLEEQGRQKVSFSFSLTKKTLQNRFLTALGNEKQSDTPNPPAVPLQVDSTPKMKMEIGDTLSTAEESSPPKSRVELGKIHFKKHLLHVTSRPLLATTTAVASPPTHAAPLPAVIAESTTVDSPPSSPPPPPPPAQATTLSSPAPVTEPVALPHTPITVLMAAPVPLPVDVAVRSLKEPPIIIVPESLEADTKQDTISNSLEEHVTQILNEQADISSKKEDSHIGKDEEIPDSSKISLSCKKTGSKKKSSQSEGIFLGSESDEDSVRTSSSQRSHDLKFSASIEKERDFKKSSAPLKSEDLGKPSRSKTDRDDKYFSYSKLERDTRYVSSRCRSERERRRSRSHSRSERGSRTNLSYSRSERSHYYDSDRRYHRSSPYRERTRYSRPYTDNRARESSDSEEEYKKTYSRRTSSHSSSYRDLRTSSYSKSDRDCKTETSYLEMERRGKYSSKLERESKRTSENEAIKRCCSPPNELGFRRGSSYSKHDSSASRYKSTLSKPIPKSDKFKNSFCCTELNEEIKQSHSFSLQTPCSKGSELRMINKNPEREKAGSPAPSNRLNDSPTLKKLDELPIFKSEFITHDSHDSIKELDSLSKVKNDQLRSFCPIELNINGSPGAESDLATFCTSKTDAVLMTSDDSVTGSELSPLVKACMLSSNGFQNISRCKEKDLDDTCMLHKKSESPFRETEPLVSPHQDKLMSMPVMTVDYSKTVVKEPVDTRVSCCKTKDSDIYCTLNDSNPSLCNSEAENIEPSVMKISSNSFMNVHLESKPVICDSRNLTDHSKFACEEYKQSIGSTSSASVNHFDDLYQPIGSSGIASSLQSLPPGIKVDSLTLLKCGENTSPVLDAVLKSKKSSEFLKHAGKETIVEVGSDLPDSGKGFASRENRRNNGLSGKCLQEAQEEGNSILPERRGRPEISLDERGEGGHVHTSDDSEVVFSSCDLNLTMEDSDGVTYALKCDSSGHAPEIVSTVHEDYSGSSESSNDESDSEDTDSDDSSIPRNRLQSVVVVPKNSTLPMEETSPCSSRSSQSYRHYSDHWEDERLESRRHLYEEKFESIASKACPQTDKFFLHKGTEKNPEISFTQSSRKQIDNRLPELSHPQSDGVDSTSHTDVKSDPLGHPNSEETVKAKIPSRQQEELPIYSSDFEDVPNKSWQQTTFQNRPDSRLGKTELSFSSSCEIPHVDGLHSSEELRNLGWDFSQEKPSTTYQQPDSSYGACGGHKYQQNAEQYGGTRDYWQGNGYWDPRSGRPPGTGVVYDRTQGQVPDSLTDDREEEENWDQQDGSHFSDQSDKFLLSLQKDKGSVQAPEISSNSIKDTLAVNEKKDFSKNLEKNDIKDRGPLKKRRQEIESDSESDGELQDRKKVRVEVEQGETSVPPGSALV""".replace('\n', '')
    fasta_SETD2_IDR2 = """VYRIPKKSQTEKENTTTERGRDAVGFRDQTPAPKTPNRSRERDPDKQTQNKEKRKRRSSLSPPSSAYERGTKRPDDRYDTPTSKKKVRIKDRNKLSTEERRKLFEQEVAQREAQKQQQQMQNLGMTSPLPYDSLGYNAPHHPFAGYPPGYPMQAYVDPSNPNAGKVLLPTPSMDPVCSPAPYDHAQPLVGHSTEPLSAPPPVPVVPHVAAPVEVSSSQYVAQSDGVVHQDSSVAVLPVPAPGPVQGQNYSVWDSNQQSVSVQQQYSPAQSQATIYYQGQTCPTVYGVTSPYSQTTPPIVQSYAQPSLQYIQGQQIFTAHPQGVVVQPAAAVTTIVAPGQPQPLQPSEMVVTNNLLDLPPPSPPKPKTIVLPPNWKTARDPEGKIYYYHVITRQTQWDPPTWESPGDDASLEHEAEMDLGTPTYDENPMKASKKPKTAEAD
""".replace('\n', '')
    fasta_SCAF1_IDR1 = """MEEEDESRGKTEESGEDRGDGPPDRDPTLSPSAFILRAIQQAVGSSLQGDLPNDKDGSRCHGLRWRRCRSPRSEPRSQESGGTDTATVLDMATDSFLAGLVSVLDPPDTWVPSRLDLRPGESEDMLELVAEVRIGDRDPIPLPVPSLLPRLRAWRTGKTVSPQSNSSRPTCARHLTLGTGDGGPAPPPAPSSASSSPSPSPSSSSPSPPPPPPPPAPPAPPAPRFDIYDPFHPTDEAYSPPPAPEQKYDPFEPTGSNPSSSAGTPSPEEEEEEEEEEEEEEEDEEEEEGLSQSISRISETLAGIYDDNSLSQDFPGDESPRPDAQPTQPTPAPGTPPQVDSTRADGAMRRRVFVVGTEAEACREGKVSVEVVTAGGAALPPPLLPPGDSEIEEGEIVQPEEEPRLALSLFRPGGRAARPTPAASATPTAQPLPQPPAPRAPEGDDFLSLHAESDGEGALQVDLGEPAPAPPAADSRWGGLDLRRKILTQRRERYRQRSPSPAPAPAPAAAAGPPTRKKSRRERKRSGEAKEAASSSSGTQPAPPAPASPWDSKKHRSRDRKPGSHASSSARRRSRSRSRSRSTRRRSRSTDRRRGGSRRSRSREKRRRRRRSASPPPATSSSSSSRRERHRGKHRDGGGSKKKKKRSRSRGEKRSGDGSEKAPAPAPPPSGSTSCGDRDSRRRGAVPPSIQDLTDHDLFAIKRTITVGRLDKSDPRGPSPAPASSPKREVLYDSEGLSGEERGGKSSQKDRRRSGAASSSSSSREKGSRRKALDGGDRDRDRDRDRDRDRSSKKARPPKESAPSSGPPPKPPVSSGSGSSSSSSSCSSRKVKLQSKVAVLIREGVSSTTPAKDAASAGLGSIGVKFSRDRESRSPFLKPDERAPTEMAKAAPGSTKPKKTKVKAKAGAKKTKGTKGKTKPSKTRKKVRSGGGSGGSGGQVSLKKSKADSCSQAAGTKGAEETSWSGEERAAKVPSTPPPKAAPPPPALTPDSQTVDSSCKTPEVSFLPEEATEEAGVRGGAEEEEEEEEEEEEEEEEEEQQPATTTATSTAAAAPSTAPSAGSTAGDSGAEDGPASRVSQLPTLPPPMPWNLPAGVDCTTSGVLALTALLFKMEEANLASRAKAQELIQATNQILSHRKPPSSLGMTPAPVPTSLGLPPGPSSYLLPGSLPLGGCGSTPPTPTGLAATSDKREGSSSSEGR""".replace('\n', '')
    fasta_SRC3_IDR1 = """EQNGYRPNPNPVGQGIRPPMAGCNSSVGGMSMSPNQGLQMPSSRAYGLADPSTTGQMSGARYGGSSNIASLTPGPGMQSPSSYQNNNYGLNMSSPPHGSPGLAPNQQNIMISPRNRGSPKIASHQFSPVAGVHSPMASSGNTGNHSFSSSSLSALQAISEGVGTSLLSTLSSPGPKLDNSPNMNITQPSKVSNQDSKSPLGFYCDQNPVESSMCQSNSRDHLSDKESKESSVEGAENQRGPLESKGHKKLLQLLTCSSDDRGHSSLTNSPLDSSCKESSVSVTSPSGVSSSTSGGVSSTSNMHGSLLQEKHRILHKLLQNGNSPAEVAKITAEATGKDTSSITSCGDGNVVKQEQLSPKKKENNALLRYLLDRDDPSDALSKELQPQVEGVDNKMSQCTSSTIPSSSQEKDPKIKTETSEEGSGDLDNLDAILGDLTSSDFYNNSISSNGSHLGTKQQVFQGTNSLGLKSSQSVQSIRPPYNRAVSLDSPVSVGSSPPVKNISAFPMLPKQPMLGGNPRMMDSQENYGSSMGGPNRNVTVTQTPSSGDWGLPNSKAGRMEPMNSNSMGRPGGDYNTSLPRPALGGSIPTLPLRSNSIPGARPVLQQQQQMLQMRPGEIPMGMGANPYGQAAASNQLGSWPDGMLSMEQVSHGTQNRPLLRNSLDDLVGPPSNLEGQSDERALLDQLHTLLSNTDATGLEEIDRALGIPELVNQGQALEPKQDAFQGQEAAVMMDQKAGLYGQTYPAQGPPMQGGFHLQGQSPSFNSMMNQMNQQGNFPLQGMHPRANIMRPRTNTPKQLRMQLQQRLQGQQFLNQSRQALELKMENPTAGGAAVMRPMMQPQVSSQQGFLNAQMVAQRSRELLSHHFRQQRVAMMMQQQQQQQQQQQQQQQQQQQQQQQQQQQQQTQAFSPPPNVTASPSMDGLLAGPTMPQAPPQQFPYQPNYGMGQQPDPAFGRVSSPPNAMMSSRMGPSQNPMMQHPQAASIYQSSEMKGWPSGNLARNSSFSQQQFAHQGNPAVYSMVHMNGSSGHMGQMNMNPMPMSGMPMGPDQKY""".replace('\n', '')
    fasta_Gli3_IDR1 = """MEAQSHSSTTTEKKKVENSIVKCSTRTDVSEKAVASSTTSNEDESPGQTYHRERRNAITMQPQNVQGLSKVSEEPSTSSDERASLIKKEIHGSLPHVAEPSVPYRGTVFAMDPRNGYMEPHYHPPHLFPAFHPPVPIDARHHEGRYHYDPSPIPPLHMTSALSSSPTYPDLPFIRISPHRNPTAASESPFSPPHPYINPYMDYIRSLHSSPSLSMISATRGLSPTDAPHAGVSPAEYYHQMALLTGQRSPYADIIPSAATAGTGAIHMEYLHAMDSTRFSSPRLSARPSRKRTLSISPLSDHSFDLQTMIRTSPNSLVTILNNSRSSSSASGSYGHLSASAISPALSFTYSSAPVSLHMHQQILSRQQSLGSAFGHSPPLIHPAPTFPTQRPIPGIPTVLNPVQVSSGPSESSQNKPTSESAVSSTGDPMHNKRSKIKPDEDLPSPGARGQQEQPEGTTLVKEEGDKDESKQEPEVI""".replace('\n', '')
    fasta_Gli3_IDR2 = """PEAHVTKKQRGDIHPRPPPPRDSGSHSQSRSPGRPTQGALGEQQDLSNTTSKREECLQVKTVKAEKPMTSQPSPGGQSSCSSQQSPISNYSNSGLELPLTDGGSIGDLSAIDETPIMDSTISTATTALALQARRNPAGTKWMEHVKLERLKQVNGMFPRLNPILPPKAPAVSPLIGNGTQSNNTCSLGGPMTLLPGRSDLSGVDVTMLNMLNRRDSSASTISSAYLSSRRSSGISPCFSSRRSSEASQAEGRPQNVSVADSYDPISTDASRRSSEASQSDGLPSLLSLTPAQQYRLKAKYAAATGGPPPTPLPNMERMSLKTRLALLGDALEPGVALPPVHAPRRCSDGGAHGYGRRHLQPHDAPGHGVRRASDPVRTGSEGLALPRVPRFSSLSSCNPPAMATSAEKRSLVLQNYTRPEGGQSRNFHSSPCPPSITENVTLESLTMDADANLNDEDFLPDDVVQYLNSQNQAGYEQHFPSALPDDSKVPHGPGDFDAPGLPDSHAGQQFHALEQPCPEGSKTDLPIQWNEVSSGSADLSSSKLKCGPRPAVPQTRAFGFCNGMVVHPQNPLRSGPAGGYQTLGENSNPYGGPEHLMLHNSPGSGTSGNAFHEQPCKAPQYGNCLNRQPVAPGALDGACGAGIQASKLKSTPMQGSGGQLNFGLPVAPNESAGSMVNGMQNQDPVGQGYLAHQLLGDSMQHPGAGRPGQQMLGQISATSHINIYQGPESCLPGAHGMGSQPSSLAVVRGYQPCASFGGSRRQAMPRDSLALQSGQLSDTSQTCRVNGIKMEMKGQPHPLCSNLQNYSGQFYDQTVGFSQQDTKAGSFSISDASCLLQGTSAKNSELLSPGANQVTSTVDSLDSHDLEGVQIDFDAIIDDGDHSSLMSGALSPSIIQNLSHSSSRLTTPRASLPFPALSMSTTNMAIGDMSSLLTSLAEESKFLAVM""".replace('\n', '')
    fasta_Gli2_IDR2 = """PDAHVTKKQRNDVHLRTPLLKENGDSEAGTEPGGPESTEASSTSQAVEDCLHVRAIKTESSGLCQSSPGAQSSCSSEPSPLGSAPNNDSGVEMPGTGPGSLGDLTALDDTPPGADTSALAAPSAGGLQLRKHMTTMHRFEQLKKEKLKSLKDSCSWAGPTPHTRNTKLPPLPGSGSILENFSGSGGGGPAGLLPNPRLSELSASEVTMLSQLQERRDSSTSTVSSAYTVSRRSSGISPYFSSRRSSEASPLGAGRPHNASSADSYDPISTDASRRSSEASQCSGGSGLLNLTPAQQYSLRAKYAAATGGPPPTPLPGLERMSLRTRLALLDAPERTLPAGCPRPLGPRRGSDGPTYGHGHAGAAPAFPHEAPGGGARRASDPVRRPDALSLPRVQRFHSTHNVNPGPLPPCADRRGLRLQSHPSTDGGLARGAYSPRPPSISENVAMEAVAAGVDGAGPEADLGLPEDDLVLPDDVVQYIKAHASGALDEGTGQVYPTESTGFSDNPRLPSPGLHGQRRMVAADSNVGPSAPMLGGCQLGFGAPSSLNKNNMPVQWNEVSSGTVDALASQVKPPPFPQGNLAVVQQKPAFGQYPGYSPQGLQASPGGLDSTQPHLQPRSGAPSQGIPRVNYMQQLRQPVAGSQCPGMTTTMSPHACYGQVHPQLSPSTISGALNQFPQSCSNMPAKPGHLGHPQQTEVAPDPTTMGNRHRELGVPDSALAGVPPPHPVQSYPQQSHHLAASMSQEGYHQVPSLLPARQPGFMEPQTGPMGVATAGFGLVQPRPPLEPSPTGRHRGVRAVQQQLAYARATGHAMAAMPSSQETAEAVPKGAMGNMGSVPPQPPPQDAGGAPDHSMLYYYGQIHMYEQDGGLENLGSCQVMRSQPPQPQACQDSIQPQPLPSPGVNQVSSTVDSQLLEAPQIDFDAIMDDGDHSSLFSGALSPSLLHSLSQNSSRLTTPRNSLTLPSIPAGISNMAVGDMSSMLTSLAEESKFLNMM""".replace('\n', '')
    
    proteins.loc['SETD2_IDR1'] = dict(eps_factor=0.2,pH=7.2,fasta=list(fasta_SETD2_IDR1),ionic=0.15)
    proteins.loc['SETD2_IDR2'] = dict(eps_factor=0.2,pH=7.2,fasta=list(fasta_SETD2_IDR2),ionic=0.15)
    proteins.loc['SCAF1_IDR1'] = dict(eps_factor=0.2,pH=7.2,fasta=list(fasta_SCAF1_IDR1),ionic=0.15)
    proteins.loc['SRC3_IDR1'] = dict(eps_factor=0.2,pH=7.2,fasta=list(fasta_SRC3_IDR1),ionic=0.15)
    proteins.loc['Gli3_IDR1'] = dict(eps_factor=0.2,pH=7.2,fasta=list(fasta_Gli3_IDR1),ionic=0.15)
    proteins.loc['Gli3_IDR2'] = dict(eps_factor=0.2,pH=7.2,fasta=list(fasta_Gli3_IDR2),ionic=0.15)
    proteins.loc['Gli2_IDR2'] = dict(eps_factor=0.2,pH=7.2,fasta=list(fasta_Gli2_IDR2),ionic=0.15)

    return proteins

def genTop(residues,fasta,path,L):
    N_res = len(fasta)
    top = md.Topology()
    chain = top.add_chain()
    for resname in fasta:
        residue = top.add_residue(residues.loc[resname,'three'], chain)
        top.add_atom(residues.loc[resname,'three'], element=md.element.carbon, residue=residue)
    for i in range(N_res-1):
        top.add_bond(top.atom(i),top.atom(i+1))
    pos = [[0,0,(i-N_res/2.)*.38] for i in range(N_res)]
    t = md.Trajectory(np.array(pos).reshape(N_res,3), top, 0, [L,L,L], [90,90,90])
    t.save_pdb(path+'/top.pdb')

def genParams(r,seq,temp,ionic):
    RT = 8.3145*temp*1e-3
    pH = 7.2
    fepsw = lambda T : 5321/T+233.76-0.9297*T+0.1417*1e-2*T*T-0.8292*1e-6*T**3
    epsw = fepsw(temp)
    lB = 1.6021766**2/(4*np.pi*8.854188*epsw)*6.022*1000/RT
    # Calculate the inverse of the Debye length
    yukawa_kappa = np.sqrt(8*np.pi*lB*ionic*6.022/10)
    fasta = list(seq.fasta)
    r.loc['H','q'] = 1. / ( 1 + 10**(pH-6) )
    r.loc['X'] = r.loc[fasta[0]]
    r.loc['X','q'] = r.loc[fasta[0],'q'] + 1.
    r.loc['X','MW'] = r.loc[fasta[0],'MW'] + 2.
    fasta[0] = 'X'
    r.loc['Z'] = r.loc[fasta[-1]]
    r.loc['Z','q'] = r.loc[fasta[-1],'q'] - 1.
    r.loc['Z','MW'] = r.loc[fasta[-1],'MW'] + 16.
    fasta[-1] = 'Z'
    # Calculate the prefactor for the Yukawa potential
    qq = pd.DataFrame(r.q.values*r.q.values.reshape(-1,1),index=r.q.index,columns=r.q.index)
    yukawa_eps = qq*lB*RT
    types = list(np.unique(fasta))
    pairs = np.array(list(itertools.combinations_with_replacement(types,2)))
    return yukawa_kappa, yukawa_eps, types, pairs, fasta, r

#Sequence descriptors
def calc_scd_shd(aa_params,seq,pH=7.2,beta=-1): #maybe pH should be 7.4 like above...
    fasta = list(seq.fasta)
    # set histidine charge based on pH; in sims HIS is neutral
    aa_params.loc['H','q'] = 1. / ( 1 + 10**(pH-6) )
    # set lambda parameters to a given model
    N = len(fasta)
    pairs = np.array(list(itertools.combinations(fasta,2)))
    pairs_indices = np.array(list(itertools.combinations(range(N),2)))
    # calculate sequence separations
    ij_dist = np.diff(pairs_indices,axis=1).flatten().astype(float)
    # calculate charge products
    qq = aa_params.q.loc[pairs[:,0]].values*aa_params.q.loc[pairs[:,1]].values
    # calculate lambda sums
    ll = aa_params.lambdas.loc[pairs[:,0]].values+aa_params.lambdas.loc[pairs[:,1]].values
    scd = np.sum(qq*np.sqrt(ij_dist))/N
    shd = np.sum(ll*np.power(np.abs(ij_dist),beta))/N
    ncpr = aa_params.q.loc[fasta].values.mean()
    q_sum = aa_params.q.loc[fasta].values.sum()
    #return scd,shd,ncpr,q_sum
    df_descrip = pd.DataFrame(index=['ncpr','scd','shd','q_sum','ij'],columns=['value'])
    df_descrip.loc['ncpr','value'] = ncpr
    df_descrip.loc['scd','value'] = scd
    df_descrip.loc['shd','value'] = shd
    df_descrip.loc['q_sum','value'] = q_sum
    df_descrip.loc['ij','value'] = ij_dist
    return df_descrip

#energy maps
HALR = lambda r,s,l : 4*0.8368*l*((s/r)**12-(s/r)**6)
HASR = lambda r,s,l : 4*0.8368*((s/r)**12-(s/r)**6)+0.8368*(1-l)
HA = lambda r,s,l : np.where(r<2**(1/6)*s, HASR(r,s,l), HALR(r,s,l))
HASP = lambda r,s,l,rc : np.where(r<rc, HA(r,s,l)-HA(rc,s,l), 0)

def calcEnergyMap(t,df,prot,rc):
    indices = t.top.select_pairs('all','all')
    mask = np.abs(indices[:,0]-indices[:,1])>1 #exclude >1, was used to exclude bonded pairs
    indices = indices[mask]
    d = md.compute_distances(t,indices) #distances between pairs for each frame
    # d[d>rc] = np.inf
    pairs = np.array(list(combinations(list(prot.fasta),2)))
    pairs = pairs[mask]
    sigmas = 0.5*(df.loc[pairs[:,0]].sigmas.values+df.loc[pairs[:,1]].sigmas.values)
    lambdas = 0.5*(df.loc[pairs[:,0]].lambdas.values+df.loc[pairs[:,1]].lambdas.values)
    emap = np.zeros(pairs.shape[0])
    print(d.shape,d.shape[0]//20*20)
    d = d[:d.shape[0]//20*20]
    for i,r in enumerate(np.split(d,20,axis=0)):
        emap += np.nansum(HASP(r,sigmas[np.newaxis,:],lambdas[np.newaxis,:],rc),axis=0)
    return indices, emap/d.shape[0]

def calcRg(t,residues,seq):
    fasta = list(seq.fasta)
    masses = residues.loc[fasta,'MW'].values
    # calculate the center of mass
    cm = np.sum(t.xyz*masses[np.newaxis,:,np.newaxis],axis=1)/masses.sum()
    # calculate residue-cm distances
    si = np.linalg.norm(t.xyz - cm[:,np.newaxis,:],axis=2)
    # calculate rg
    rgarray = np.sqrt(np.sum(si**2*masses,axis=1)/masses.sum())
    return rgarray

def calcRs(traj):
    pairs = traj.top.select_pairs('all','all')
    d = md.compute_distances(traj,pairs)
    dmean = d.mean(axis=0)
    nres = traj.n_atoms
    ij = np.arange(2,nres,1)
    diff = [x[1]-x[0] for x in pairs]
    dij = np.empty(0)
    for i in ij:
        dij = np.append(dij,dmean[diff==i].mean())
    ln_ij =  np.log(ij[ij>10])
    ln_dij = np.log(dij[ij>10])
    return ij,dij,ln_ij,ln_dij,np.mean(1/d,axis=1)

def analyse(residues,path,seq):
    top = md.Topology()
    chain = top.add_chain()
    if os.path.exists(path+'/traj.gsd'):
        traj = md.load(path+'/traj.gsd')
    else:
        traj = md.load_xtc(path+'/traj.xtc',top=path+'/top.pdb')
    N_res = traj.n_atoms
    fasta = list(seq.fasta)
    #fixing the trajectory to the middle of the box
    for resname in fasta:
        residue = top.add_residue(residues.loc[resname,'three'],chain)
        top.add_atom(residues.loc[resname,'three'], element=md.element.carbon, residue=residue)
    for i in range(N_res-1):
        top.add_bond(top.atom(i),top.atom(i+1))
    traj.top = top
    traj = traj.image_molecules(inplace=False, anchor_molecules=[set(traj.top.atoms)],make_whole=True)
    print('Number of frames: {:d}'.format(traj.n_frames))
    if os.path.exists(path+'/traj.gsd'):
        traj[-1].save_pdb(path+'/top.pdb')
        traj.save_xtc(path+'/traj.xtc')
        os.remove(path+'/traj.gsd')
    #skip first 10 frames
    traj = traj[10:]
    #energy maps
    df_map = pd.DataFrame(index=range(traj.n_atoms),columns=range(traj.n_atoms),dtype=float)
    pairs, emap = calcEnergyMap(traj,residues,seq,2.0)
    for k,(i,j) in enumerate(pairs):
        df_map.loc[i,j] = emap[k]
        df_map.loc[j,i] = emap[k]
    df_analysis = pd.DataFrame(index=['Rg','ete','rh','nu','R0','nu_ln','R0_ln'],columns=['value','error'])
    #rg
    rgarray = calcRg(traj,residues,seq)
    np.save(path+'/rg.npy',rgarray)
    block_rg = BlockAnalysis(rgarray, multi=1)
    block_rg.SEM()
    df_analysis.loc['Rg','value'] = block_rg.av
    df_analysis.loc['Rg','error'] = block_rg.sem
    #ete
    ete = md.compute_distances(traj,atom_pairs=[[0,N_res-1]]).flatten()
    np.save(path+'/ete.npy',ete)
    block_ete = BlockAnalysis(ete, multi=1)
    block_ete.SEM()
    df_analysis.loc['ete','value'] = block_ete.av
    df_analysis.loc['ete','error'] = block_ete.sem
    #nonlinear scaling exponent
    f = lambda x,R0,v : R0*np.power(x,v)
    ij,dij,ln_ij,ln_dij,invrij = calcRs(traj)
    block_invrij = BlockAnalysis(invrij, multi=1)
    block_invrij.SEM()
    df_analysis.loc['rh','value'] = 1/(1-1/N_res)/block_invrij.av
    df_analysis.loc['rh','error'] = block_invrij.sem/(1-1/N_res)/block_invrij.av/block_invrij.av
    np.save(path+'/rs.npy',dij)
    popt, pcov = curve_fit(f,ij[ij>10],dij[ij>10],p0=[.4,.5])
    df_analysis.loc['nu','value'] = popt[1]
    df_analysis.loc['nu','error'] = pcov[1,1]**0.5
    df_analysis.loc['R0','value'] = popt[0]
    df_analysis.loc['R0','error'] = pcov[0,0]**0.5
    #linear scaling exponent
    f = lambda x,R0,v : R0 + v*x
    popt, pcov = curve_fit(f, ln_ij, ln_dij, p0=[-1,.4],bounds=([-3,0],[1,1]))
    df_analysis.loc['nu_ln','value'] = popt[1]
    df_analysis.loc['nu_ln','error'] = pcov[1,1]**0.5
    df_analysis.loc['R0_ln','value'] = np.exp(popt[0])
    df_analysis.loc['R0_ln','error'] = np.exp(popt[0])*(pcov[0,0]**0.5) #error propagation
    return df_map,df_analysis

def xy_spiral_array(n, delta=0, arc=.38, separation=.7):
    """
    create points on an Archimedes' spiral
    with `arc` giving the length of arc between two points
    and `separation` giving the distance between consecutive 
    turnings
    """
    def p2c(r, phi):
        """
        polar to cartesian
        """
        return (r * np.cos(phi), r * np.sin(phi))
    r = arc
    b = separation / (2 * np.pi)
    phi = float(r) / b
    coords = []
    for i in range(n):
        coords.append(list(p2c(r, phi))+[0])
        phi += float(arc) / r
        r = b * phi
    return np.array(coords)+delta

def simulate(residues,sequences,seq_name,path):
    hoomd.context.initialize("--mode=gpu");
    hoomd.option.set_notice_level(1)
    hoomd.util.quiet_status()

    seq = sequences.loc[seq_name]
    df_descrip = calc_scd_shd(residues,seq)
    df_descrip.to_pickle(path+'/descriptors.pkl')
    df_descrip.to_csv(path+'/descriptors.csv')

    lj_eps = 4.184*.2
    temp = 298
    ionic_strength = 0.15 # M
    RT = 8.3145*temp*1e-3

    yukawa_kappa, yukawa_eps, types, pairs, fasta, residues = genParams(residues,seq,temp,ionic_strength)

    sigmamap = pd.DataFrame((residues.sigmas.values+residues.sigmas.values.reshape(-1,1))/2,
                            index=residues.sigmas.index,columns=residues.sigmas.index)
    lambdamap = pd.DataFrame((residues.lambdas.values+residues.lambdas.values.reshape(-1,1))/2,
                            index=residues.lambdas.index,columns=residues.lambdas.index)

    N_res = seq.N
    L = 200 #(N_res-1)*0.38+4
    N_save = 7000 if N_res < 150 else int(np.ceil(3e-4*N_res**2)*1000)
    N_steps = 1010*N_save

    genTop(residues,fasta,path,L)

    snapshot = hoomd.data.make_snapshot(N=N_res,
                                box=hoomd.data.boxdim(Lx=L, Ly=L, Lz=L),
                                particle_types=types,
                                bond_types=['polymer']);

    snapshot.bonds.resize(N_res-1);

    snapshot.particles.position[:] = xy_spiral_array(N_res) #[[0,0,(i-N_res/2.)*.38] for i in range(N_res)]
    snapshot.particles.typeid[:] = [types.index(a) for a in fasta]
    snapshot.particles.mass[:] = [residues.loc[a].MW for a in fasta]

    snapshot.bonds.group[:] = [[i,i+1] for i in range(N_res-1)];
    snapshot.bonds.typeid[:] = [0] * (N_res-1)

    hoomd.init.read_snapshot(snapshot);

    hb = hoomd.md.bond.harmonic();
    hb.bond_coeff.set('polymer', k=8033.0, r0=0.38);

    nl = hoomd.md.nlist.cell();

    ah = azplugins.pair.ashbaugh(r_cut=2.0, nlist=nl)
    yukawa = hoomd.md.pair.yukawa(r_cut=4.0, nlist=nl)
    for a,b in pairs:
        ah.pair_coeff.set(a, b, lam=lambdamap.loc[a,b], epsilon=lj_eps, sigma=sigmamap.loc[a,b], r_cut=2.0)
        yukawa.pair_coeff.set(a, b, epsilon=yukawa_eps.loc[a,b], kappa=yukawa_kappa, r_cut=4.)

    ah.set_params(mode='shift')
    yukawa.set_params(mode='shift')
    nl.reset_exclusions(exclusions = ['bond'])

    integrator_mode = hoomd.md.integrate.mode_standard(dt=0.005);
    integrator = hoomd.md.integrate.langevin(group=hoomd.group.all(),kT=RT,seed=np.random.randint(100));

    for a in types:
        integrator.set_gamma(a, residues.loc[a].MW/100)

    hoomd.run(10000)
    integrator.disable()

    integrator_mode = hoomd.md.integrate.mode_standard(dt=0.01);
    integrator = hoomd.md.integrate.langevin(group=hoomd.group.all(),kT=RT,seed=np.random.randint(100));

    for a in types:
        integrator.set_gamma(a, residues.loc[a].MW/100)

    hoomd.dump.gsd(filename=path+'/traj.gsd', period=N_save, group=hoomd.group.all(), overwrite=True);

    hoomd.run(N_steps)

    hoomd.dump.gsd(filename=path+'/restart.gsd', group=hoomd.group.all(), truncate=True, period=None, phase=0)

residues = pd.read_csv('residues.csv').set_index('one',drop=False)

sequences = initProteins()
sequences['N'] = sequences['fasta'].apply(lambda x : len(x))
sequences.to_pickle('proteins.pkl')

t0 = time.time()
simulate(residues,sequences,args.seq_name,args.path)
df_map,df_analysis = analyse(residues,args.path,sequences.loc[args.seq_name])
df_analysis.to_csv(args.path+'/analysis.csv')
df_map.to_csv(args.path+'/map.csv')
print('Timing sim and analysis {:.3f}'.format((time.time()-t0)/3600))
