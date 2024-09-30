#include "cut_right_lobe.h"
#include "liver_cut.h"  
#include <vtkSTLWriter.h>
#include <iostream>

// 定义函数，执行第二刀的操作
std::tuple<vtkSmartPointer<vtkPolyData>,vtkSmartPointer<vtkPolyData>> cut_right_lobe(vtkSmartPointer<vtkPolyData> liver, 
                    vtkSmartPointer<vtkPolyData> rb1v, 
                    vtkSmartPointer<vtkPolyData> rb2v, 
                    vtkSmartPointer<vtkPolyData> rb3v, 
                    vtkSmartPointer<vtkPolyData> rf_up_v, 
                    vtkSmartPointer<vtkPolyData> rf_down_v, 
                    vtkSmartPointer<vtkPolyData> left_inside_v, 
                    vtkSmartPointer<vtkPolyData> tail,
                    double (&mbounds)[6]) 
{
    // 第2刀，提取右半叶
    auto exc_left_inside_v = mergePolyData(rf_up_v, rf_down_v, rb1v, rb2v, rb3v);
    auto left_side = mergePolyData(left_inside_v, tail);
    auto [cut_left_inside, cut_left_inside_points] = find_points(exc_left_inside_v, left_side);
    auto imp2 = plane_generator(cut_left_inside);
    auto [liver_right, liver_left_pre] = cut_run2(imp2, liver,mbounds);

    // 保存右半叶切割后的 STL 文件
    vtkSmartPointer<vtkSTLWriter> stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
    stlWriter->SetFileName("C:/code/liver_cut/new_stl/liver_right.stl");
    stlWriter->SetInputData(liver_right);
    stlWriter->Write();

    vtkSmartPointer<vtkSampleFunction> sampleFunction = vtkSmartPointer<vtkSampleFunction>::New();
    sampleFunction->SetImplicitFunction(imp2); // 您的隐式函数
    sampleFunction->SetModelBounds(mbounds[0], mbounds[1], mbounds[2], mbounds[3], mbounds[4], mbounds[5]); // 设置采样范围
    sampleFunction->SetSampleDimensions(50, 50, 50); // 设置采样分辨率
    sampleFunction->ComputeNormalsOff();
    sampleFunction->Update();

    vtkSmartPointer<vtkContourFilter> contourFilter = vtkSmartPointer<vtkContourFilter>::New();
    contourFilter->SetInputConnection(sampleFunction->GetOutputPort());
    contourFilter->SetValue(0, 0.0); // 提取等值为 0 的面，即隐式函数的零等值面
    contourFilter->Update();

    stlWriter->SetFileName("C:/code/liver_cut/new_imp/imp_right_lobe.stl"); // 设置输出文件名
    stlWriter->SetInputConnection(contourFilter->GetOutputPort());
    stlWriter->Write();


    std::cout << "Successfully saved liver_right cut" << std::endl;
    std::cout << "Liver_right cut completed" << std::endl;

    return std::make_tuple(liver_right,liver_left_pre);
}
