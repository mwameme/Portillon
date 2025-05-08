#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <random>
#include <complex>
#include <cmath> 



using namespace std;
using namespace Eigen;



int main()
{
    double lambda, mu;
    int N;
    double temps_arret;

    cout << "Valeur de lambda (portillon)" << endl;
    cin >> lambda;
    cout << "Valeur de mu" << endl;
    cin >> mu;
    cout << "taille matrice (> 3)" << endl;
    cin >> N;
    if (N < 3)
        N = 3;
    cout << "temps de simulation" << endl;
    cin >> temps_arret;

    // vector<vector<double>> matrice_iter(N, vector<double>(N, 0.));
    // matrice_iter[0][1] = 1.;
    // matrice_iter[1][0] = lambda / (lambda + mu);
    // matrice_iter[1][2] = mu / (lambda + mu);
    // matrice_iter[2][1] = 2.0 * lambda / (2.0 * lambda + mu);
    // matrice_iter[2][3] = mu / (2.0 * lambda + mu);

    // for (int i(3); i < N - 1; ++i) {
    //     matrice_iter[i][i - 1] = 3.0 * lambda / (3.0 * lambda + mu);
    //     matrice_iter[i][i + 1] = mu / (3.0 * lambda + mu);
    // }
    // matrice_iter[N - 1][N - 2] = 3.0 * lambda / (3.0 * lambda + mu);
    // matrice_iter[N - 1][N - 1] = 1 - ( 3.0 * lambda / (3.0 * lambda + mu) );

    vector<vector<double>> matrice_continue(N, vector<double>(N, 0.));
    matrice_continue[0][0] = -mu;
    matrice_continue[0][1] = mu;
    matrice_continue[1][0] = lambda;
    matrice_continue[1][1] = -lambda - mu;
    matrice_continue[1][2] = mu;
    matrice_continue[2][1] = 2.0 * lambda;
    matrice_continue[2][2] = -2.0 * lambda - mu;
    matrice_continue[2][3] = mu;

    for (int i(3); i < (N - 1); ++i) {
        matrice_continue[i][i - 1] = 3.0 * lambda;
        matrice_continue[i][i + 1] = mu;
        matrice_continue[i][i] = -3.0 * lambda - mu;
    }
    matrice_continue[N - 1][N - 2] = 3.0 * lambda;
    matrice_continue[N - 1][N - 1] = -3.0 * lambda;


    vector<double> temps_liste(0);
    vector<double> taille_liste(0);
    random_device rd;                     // Source d'entropie
    mt19937 gen(rd());                    // Générateur Mersenne Twister
    exponential_distribution<> dist_0(mu);
    exponential_distribution<> dist_1(mu + lambda);
    exponential_distribution<> dist_2(mu + 2.0 * lambda);
    exponential_distribution<> dist_3(mu + 3.0 * lambda);
    uniform_real_distribution<> uniform(0.0, 1.0);

    double temps_local = 0.0;
    int taille = 0;
    while (temps_local < temps_arret) {
        if (taille == 0) {
            double delta_t = dist_0(gen);
            temps_local += delta_t;
            taille = 1;
        }
        else if (taille == 1) {
            double delta_t = dist_1(gen);
            double p = lambda / (lambda + mu);
            double x = uniform(gen);
            if (x < p)
                taille = 0;
            else
                taille = 2;
            temps_local += delta_t;
        }
        else if (taille == 2) {
            double delta_t = dist_2(gen);
            double p = 2.0 * lambda / (2.0 * lambda + mu);
            double x = uniform(gen);
            if (x < p)
                taille = 1;
            else
                taille = 3;
            temps_local += delta_t;
        }
        if (taille >= 3) {
            double delta_t = dist_3(gen);
            double p = 3.0 * lambda / (3.0 * lambda + mu);
            double x = uniform(gen);
            if (x < p)
                taille -= 1;
            else
                taille += 1;
            temps_local += delta_t;
        }

        temps_liste.push_back(temps_local);
        taille_liste.push_back(taille);
    }

    MatrixXcd mat(N, N);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            mat(j, i) = std::complex<double>(matrice_continue[i][j], 0.0); // transpose !
        }
    }

    ComplexEigenSolver<MatrixXcd> solver(mat);
    VectorXcd eigenvalues = solver.eigenvalues();

    double max = eigenvalues[0].real();
    int index_max = 0;
    for (int i = 1; i < eigenvalues.size(); ++i) {
        if (eigenvalues[i].real() > max) {
            max = eigenvalues[i].real();
            index_max = i;
        }
    }

    double max_2;
    int index_max_2;
    if (index_max == 0) {
        max_2 = eigenvalues[1].real();
        index_max_2 = 1;
    }
    else {
        max_2 = eigenvalues[0].real();
        index_max_2 = 0;
    }

    for (int i(0); i < eigenvalues.size(); ++i) {
        if (i == index_max) continue;
        if (eigenvalues[i].real() > max_2) {
            max_2 = eigenvalues[i].real();
            index_max_2 = i;
        }
    }

    cout << "GAP : exp(lambda_2) = " << exp(max_2) << endl;
    cout << "sup(vp) 0? = " << max << endl;

    cout << "Proba stationnaire" << endl << endl;
    VectorXcd eigenvector = solver.eigenvectors().col(index_max);
    complex<double> somme(0., 0.);
    for (int i = 0; i < eigenvector.size(); ++i)
        somme += eigenvector[i];
    for (int i = 0; i < eigenvector.size(); ++i) {
        eigenvector[i] = eigenvector[i] / somme;
        cout << eigenvector[i].real() << endl;
    }

    cout << "somme abs de la partie complexe : ";
    double somme_i = 0.;
    for (int i = 0; i < eigenvector.size(); ++i) {
        somme_i += abs(eigenvector[i].imag());
    }
    cout << somme_i << endl;

    cout << "Fin de proba stationnaire" << endl << endl;
    cout << "Simulation. Affiché : (début du temps ; taille de la file à partir de ce moment)" << endl << endl;
    for (int i = 0; i < temps_liste.size(); ++i) {
        cout << temps_liste[i] << " ; " << taille_liste[i] << endl;
    }

    cout << "Fin de simulation" << endl << endl;

}

