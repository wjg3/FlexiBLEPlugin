#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
using namespace std;

int main()
{
    fstream fout("PosVec.h", ios::out);
    fout << "#include <vector>" << endl;
    fout << endl;

    fstream fin("NAcoor.txt", ios::in);
    vector<vector<double>> NAcoor;
    for (int i = 0; i < 200; i++)
    {
        double x, y, z;
        fin >> x >> y >> z;
        vector<double> temp = {x, y, z};
        NAcoor.emplace_back(temp);
    }
    fout << "std::vector<std::vector<double>> NAPositions = {";
    for (int i = 0; i < NAcoor.size(); i++)
    {
        fout << "{" << setprecision(6) << NAcoor[i][0] << "," << NAcoor[i][1] << "," << NAcoor[i][2] << "}";
        if (i != 199)
            fout << ",";
    }
    fout << "};" << endl;

    fstream fin1("NAvel.txt", ios::in);
    vector<vector<double>> NAvel;
    for (int i = 0; i < 200; i++)
    {
        double x, y, z;
        fin1 >> x >> y >> z;
        vector<double> temp = {x, y, z};
        NAvel.emplace_back(temp);
    }
    fout << "std::vector<std::vector<double>> NAVelocities = {";
    for (int i = 0; i < NAvel.size(); i++)
    {
        fout << "{" << setprecision(6) << NAvel[i][0] << "," << NAvel[i][1] << "," << NAvel[i][2] << "}";
        if (i != 199)
            fout << ",";
    }
    fout << "};" << endl;

    fstream fin2("Neoncoor.txt", ios::in);
    vector<vector<double>> Neoncoor;
    for (int i = 0; i < 200; i++)
    {
        double x, y, z;
        fin2 >> x >> y >> z;
        vector<double> temp = {x, y, z};
        Neoncoor.emplace_back(temp);
    }
    fout << "std::vector<std::vector<double>> NeonPositions = {";
    for (int i = 0; i < Neoncoor.size(); i++)
    {
        fout << "{" << setprecision(6) << Neoncoor[i][0] << "," << Neoncoor[i][1] << "," << Neoncoor[i][2] << "}";
        if (i != 199)
            fout << ",";
    }
    fout << "};" << endl;

    fstream fin3("Neonvel.txt", ios::in);
    vector<vector<double>> Neonvel;
    for (int i = 0; i < 200; i++)
    {
        double x, y, z;
        fin3 >> x >> y >> z;
        vector<double> temp = {x, y, z};
        Neonvel.emplace_back(temp);
    }
    fout << "std::vector<std::vector<double>> NeonVelocities = {";
    for (int i = 0; i < Neonvel.size(); i++)
    {
        fout << "{" << setprecision(6) << Neonvel[i][0] << "," << Neonvel[i][1] << "," << Neonvel[i][2] << "}";
        if (i != 199)
            fout << ",";
    }
    fout << "};" << endl;

    fstream fin4("water_coor.txt", ios::in);
    vector<vector<double>> wat_coor;
    for (int i = 0; i < 600; i++)
    {
        double x, y, z;
        fin4 >> x >> y >> z;
        vector<double> temp = {x, y, z};
        wat_coor.emplace_back(temp);
    }
    fout << "std::vector<std::vector<double>> WaterPositions = {";
    for (int i = 0; i < wat_coor.size(); i++)
    {
        fout << "{" << setprecision(6) << wat_coor[i][0] << "," << wat_coor[i][1] << "," << wat_coor[i][2] << "}";
        if (i != 599)
            fout << ",";
    }
    fout << "};" << endl;

    fstream fin5("water_vel.txt", ios::in);
    vector<vector<double>> wat_vel;
    for (int i = 0; i < 600; i++)
    {
        double x, y, z;
        fin5 >> x >> y >> z;
        vector<double> temp = {x, y, z};
        wat_vel.emplace_back(temp);
    }
    fout << "std::vector<std::vector<double>> WaterVelocities = {";
    for (int i = 0; i < wat_vel.size(); i++)
    {
        fout << "{" << setprecision(6) << wat_vel[i][0] << "," << wat_vel[i][1] << "," << wat_vel[i][2] << "}";
        if (i != 599)
            fout << ",";
    }
    fout << "};" << endl;
}