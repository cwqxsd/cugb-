#pragma once

#include <Eigen/Core>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <cmath>
#include <array>
#include <chrono>
#include <vector>
#include <string>
#include <memory>

#include <rerun.hpp>
#include <behaviortree_cpp/behavior_tree.h>
#include <behaviortree_cpp/bt_factory.h>

#include "types.h"

using namespace std;
namespace chr = std::chrono;

class Locator
{
public:
    // -------------------- params --------------------
    double convergeTolerance = 0.2;
    double residualTolerance = 0.4;
    int    maxIteration = 20;
    double muOffset = 2.0;
    double numShrinkRatio = 0.85;
    double offsetShrinkRatio = 0.8;
    int    minMarkerCnt = 3;

    bool   enableLog = false;
    string logIP = "127.0.0.1:9876";

    // -------------------- state --------------------
    // NOTE: 原代码是 const RecordingStream，但 init() 里 connect/save 可能需要非 const
    rerun::RecordingStream log = rerun::RecordingStream("locator", "locator");

    vector<FieldMarker> fieldMarkers;
    // 按 type 分组后的索引：minDist/getOffset 只扫同类点
    array<vector<FieldMarker>, 256> fieldMarkersByType;

    FieldDimensions fieldDimensions;

    // hypos: [x, y, theta, residual, prob, cdf]
    using HypoMat = Eigen::Matrix<double, Eigen::Dynamic, 6, Eigen::RowMajor>;
    HypoMat hypos;
    HypoMat hypos_buf; // 双缓冲：genParticles 用它写新粒子，然后 swap

    PoseBox2D constraints;
    double offsetX = 0.0, offsetY = 0.0, offsetTheta = 0.0;

    Pose2D bestPose{0.0, 0.0, 0.0};
    double bestResidual = std::numeric_limits<double>::infinity();

    // 预计算观测 marker 的权重，避免每粒子重复 sqrt/div
    struct ObsMarker
    {
        char type;
        double x;
        double y;
        double inv_dist_w; // 3.0 / max(norm(x,y), 0.1)
    };

    void init(const FieldDimensions& fd,
              int minMarkerCnt = 4,
              double residualTolerance = 0.4,
              double muOffsetParam = 2.0,
              bool enableLog = false,
              string logIP = "127.0.0.1:9876");

    void calcFieldMarkers(const FieldDimensions& fd);

    LocateResult locateRobot(const vector<FieldMarker>& markers_r,
                             PoseBox2D constraints,
                             int numParticles = 200,
                             double offsetX = 2.0,
                             double offsetY = 2.0,
                             double offsetTheta = M_PI / 4);

    int genInitialParticles(int num = 200);
    int genParticles();

    // ---- by-type 扫描 ----
    double minDist(const FieldMarker& marker);
    array<double, 2> getOffset(const FieldMarker& marker);

    // 保持旧接口（BT 节点验证会调用）
    double residual(const vector<FieldMarker>& markers_r, const Pose2D& pose);

    bool isConverged();

    // 保持旧接口：内部会 buildObsMarkers 后走新 calcProbs(obs)
    int calcProbs(const vector<FieldMarker>& markers_r);

    Pose2D finalAdjust(const vector<FieldMarker>& markers_r, const Pose2D& pose);

    void logParticles();

private:
    inline void rebuildFieldMarkerIndex()
    {
        for (auto& v : fieldMarkersByType) v.clear();
        for (const auto& m : fieldMarkers)
        {
            fieldMarkersByType[static_cast<unsigned char>(m.type)].push_back(m);
        }
    }

    inline vector<ObsMarker> buildObsMarkers(const vector<FieldMarker>& markers_r) const;

    // 主路径用的 residual：输入已经预计算好的 ObsMarker（inv_dist_w）
    inline double residual(const vector<ObsMarker>& obs, const Pose2D& pose);

    // 主路径用的 calcProbs：Welford + 常数提取
    int calcProbs(const vector<ObsMarker>& obs);

    // 不构造 FieldMarker，直接用 type + x/y
    inline double minDistXY(char type, double x, double y) const;
    inline array<double, 2> getOffsetXY(char type, double x, double y) const;
};

/* ------------------------ Behavior Tree Nodes ------------------------ */

class Brain;
using namespace BT;

void RegisterLocatorNodes(BT::BehaviorTreeFactory &factory, Brain* brain);

class SelfLocate : public SyncActionNode
{
public:
    SelfLocate(const string &name, const NodeConfig &config, Brain *_brain) : SyncActionNode(name, config), brain(_brain) {}
    NodeStatus tick() override;

    static PortsList providedPorts()
    {
        return {
            InputPort<string>("mode", "enter_field", "must be one of [trust_direction, face_forward, fall_recovery]"),
            InputPort<double>("msecs_interval", 10000, "防止过于频繁地校准, 如果上一次校准距离现在小于这个时间, 则不重新校准."),
        };
    };

private:
    Brain *brain;
};

class SelfLocateEnterField : public SyncActionNode
{
public:
    SelfLocateEnterField(const string &name, const NodeConfig &config, Brain *_brain) : SyncActionNode(name, config), brain(_brain) {}
    NodeStatus tick() override;

    static PortsList providedPorts()
    {
        return {
            InputPort<double>("msecs_interval", 1000, "防止过于频繁地校准, 如果上一次校准距离现在小于这个时间, 则不重新校准."),
        };
    };

private:
    Brain *brain;
};

class SelfLocate1M : public SyncActionNode
{
public:
    SelfLocate1M(const string &name, const NodeConfig &config, Brain *_brain) : SyncActionNode(name, config), brain(_brain) {}
    NodeStatus tick() override;

    static PortsList providedPorts()
    {
        return {
            InputPort<double>("msecs_interval", 1000, "防止过于频繁地校准, 如果上一次校准距离现在小于这个时间, 则不重新校准."),
            InputPort<double>("max_dist", 2.0, "marker 距离机器人的距离小于此值时, 才进行校准. (距离小测距更准)"),
            InputPort<double>("max_drift", 1.0, "校准后的位置与原位置距离应小于此值, 否则认为校准失败"),
            InputPort<bool>("validate", true, "校准后, 用其它的 marker 进行验证, 要求小于 locator 的 max residual"),
        };
    };

private:
    Brain *brain;
};

class SelfLocate2X : public SyncActionNode
{
public:
    SelfLocate2X(const string &name, const NodeConfig &config, Brain *_brain) : SyncActionNode(name, config), brain(_brain) {}
    NodeStatus tick() override;

    static PortsList providedPorts()
    {
        return {
            InputPort<double>("msecs_interval", 1000, "防止过于频繁地校准, 如果上一次校准距离现在小于这个时间, 则不重新校准."),
            InputPort<double>("max_dist", 2.0, "penalty point 距离机器人的距离小于此值时, 才进行校准. (距离小测距更准)"),
            InputPort<double>("max_drift", 1.0, "校准后的位置与原位置距离应小于此值, 否则认为校准失败"),
            InputPort<bool>("validate", true, "校准后, 用其它的 marker 进行验证, 要求小于 locator 的 max residual"),
        };
    };

private:
    Brain *brain;
};

class SelfLocate2T : public SyncActionNode
{
public:
    SelfLocate2T(const string &name, const NodeConfig &config, Brain *_brain) : SyncActionNode(name, config), brain(_brain) {}
    NodeStatus tick() override;

    static PortsList providedPorts()
    {
        return {
            InputPort<double>("msecs_interval", 1000, "防止过于频繁地校准, 如果上一次校准距离现在小于这个时间, 则不重新校准."),
            InputPort<double>("max_dist", 2.0, "两个 TCross 距离机器人的距离小于此值时, 才进行校准. (距离小测距更准)"),
            InputPort<double>("max_drift", 1.0, "校准后的位置与原位置距离应小于此值, 否则认为校准失败"),
            InputPort<bool>("validate", true, "校准后, 用其它的 marker 进行验证, 要求小于 locator 的 max residual"),
        };
    };

private:
    Brain *brain;
};

class SelfLocateLT : public SyncActionNode
{
public:
    SelfLocateLT(const string &name, const NodeConfig &config, Brain *_brain) : SyncActionNode(name, config), brain(_brain) {}
    NodeStatus tick() override;

    static PortsList providedPorts()
    {
        return {
            InputPort<double>("msecs_interval", 1000, "防止过于频繁地校准, 如果上一次校准距离现在小于这个时间, 则不重新校准."),
            InputPort<double>("max_dist", 2.0, "penalty point 距离机器人的距离小于此值时, 才进行校准. (距离小测距更准)"),
            InputPort<double>("max_drift", 1.0, "校准后的位置与原位置距离应小于此值, 否则认为校准失败"),
            InputPort<bool>("validate", true, "校准后, 用其它的 marker 进行验证, 要求小于 locator 的 max residual"),
        };
    };

private:
    Brain *brain;
};

class SelfLocatePT : public SyncActionNode
{
public:
    SelfLocatePT(const string &name, const NodeConfig &config, Brain *_brain) : SyncActionNode(name, config), brain(_brain) {}
    NodeStatus tick() override;

    static PortsList providedPorts()
    {
        return {
            InputPort<double>("msecs_interval", 1000, "防止过于频繁地校准, 如果上一次校准距离现在小于这个时间, 则不重新校准."),
            InputPort<double>("max_dist", 2.0, "penalty point 距离机器人的距离小于此值时, 才进行校准. (距离小测距更准)"),
            InputPort<double>("max_drift", 1.0, "校准后的位置与原位置距离应小于此值, 否则认为校准失败"),
            InputPort<bool>("validate", true, "校准后, 用其它的 marker 进行验证, 要求小于 locator 的 max residual"),
        };
    };

private:
    Brain *brain;
};

class SelfLocateBorder : public SyncActionNode
{
public:
    SelfLocateBorder(const string &name, const NodeConfig &config, Brain *_brain) : SyncActionNode(name, config), brain(_brain) {}
    NodeStatus tick() override;

    static PortsList providedPorts()
    {
        return {
            InputPort<double>("msecs_interval", 1000, "防止过于频繁地校准, 如果上一次校准距离现在小于这个时间, 则不重新校准."),
            InputPort<double>("max_dist", 2.0, "border 距离机器人的距离小于此值时, 才进行校准. (距离小测距更准)"),
            InputPort<double>("max_drift", 1.0, "校准后的位置与原位置距离应小于此值, 否则认为校准失败"),
            InputPort<bool>("validate", true, "校准后, 用其它的 marker 进行验证, 要求小于 locator 的 max residual"),
        };
    };

private:
    Brain *brain;
};